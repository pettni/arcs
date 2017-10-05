classdef BDDSystem < handle
  %UNTITLED Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    % encoding setting (enc = {'split', 'log'})
    enc_setting = '';
    
    % BDD variable data
    s_var_num = 0;
    a_var_num = 0;
    
    % BDD encodings
    state_encodings = {};
    action_encodings = {};
    
  end
  
  properties (Access = private)
    % BDD structure ID
    % Is identifier in list of BDD systems in C
    BDD_system_ID = 0;
    
    % Extra fields used when saving object
    s_in_inds = [];
    s_out_inds = [];
    a_inds = [];
    BDD_file_string = '';
    pg_num = 0;
  end
  
  properties (Constant)
    % Encoding settings
    split_enc = 'split';
    log_enc = 'log';
  end
  
  methods
    
    function sys = BDDSystem(state_num, action_num, enc_set)
      sys.enc_setting = enc_set;
      
      % Create a bdd TransSyst with given encoding
      % but initialize states with log encoding
      sys.s_var_num = ceil(log2(state_num));
      if state_num == 1
        sys.s_var_num = 1;
      end
      for i = 1:state_num
        enc = de2bi(i-1);
        enc = [enc, zeros(1, sys.s_var_num - length(enc))];
        sys.state_encodings{i} = enc;
      end
      sys.a_var_num = ceil(log2(action_num));
      if action_num == 1
        sys.a_var_num = 1;
      end
      for i = 1:action_num
        sys.action_encodings{i} = de2bi(i-1);
      end
      
      sys.BDD_system_ID = mexBDD('initialize', sys.s_var_num, state_num,...
                                 sys.state_encodings, sys.a_var_num, ...
                                 action_num, sys.action_encodings);
    end
    
    function delete(sys)
      mexBDD('delete', sys.BDD_system_ID);
    end
    
    function add_action(sys, new_action)
      % append action with log encoding
      new_enc = de2bi(new_action-1);
      if length(new_enc) > sys.a_var_num
        sys.a_var_num = length(new_enc);
      end
      sys.action_encodings{end+1} = new_enc;
      mexBDD('add_a', sys.BDD_system_ID, new_action, sys.a_var_num, new_enc);
    end
    
    function add_state(sys, old_state, new_state)
      old_enc = sys.state_encodings{old_state};
      if strcmp(sys.enc_setting, BDDSystem.split_enc)
        % Add zero to old state index
        sys.state_encodings{old_state} = [old_enc, 0];
        % Add one to new state index
        new_enc = [old_enc, 1];
        sys.state_encodings{new_state} = new_enc;
        % Increase number of variables, if needed
        if length(old_enc)+1 > sys.s_var_num
          sys.s_var_num = length(old_enc)+1;
        end
      elseif strcmp(sys.enc_setting, sys.log_enc)
        new_enc = de2bi(new_state-1);
        sys.state_encodings{new_state} = new_enc;
        if length(new_enc) > sys.s_var_num
          sys.s_var_num = length(new_enc);
        end
      else
        disp('Invalid encoding setting!');
        return;
      end
      mexBDD('add_s', sys.BDD_system_ID, new_state, sys.s_var_num, uint32(new_enc));
    end
  
    function add_transition(sys, s1, a, s2)
      % TODO: check for invalid states or actions
%       fprintf('MATLAB: Adding transition (%d, %d, %d)\n', [s1, a, s2]');
      mexBDD('add_trans', sys.BDD_system_ID, s1, a, s2);
    end
    
    function [state1, actions, state2] = get_trans_with_state(sys, state)
      [state1, actions, state2] = mexBDD('get_trans_with_s', sys.BDD_system_ID, state);
    end
    
    function remove_trans_with_state(sys, state)
      mexBDD('rm_trans_with_s', sys.BDD_system_ID, state);
    end
    
    % Simply sets the encoding for a state index
    % WARNING: Does not increment state counter and old state enc not
    % deleted
    % Use to change id of already or change encoding for added state
    function set_state_enc(sys, ind, enc)
      mexBDD('set_state_enc', sys.BDD_system_ID, ind, enc);
      sys.state_encodings{ind} = enc;
    end
    
    function change_state(sys, ind, new_ind)
      old_enc = sys.state_encodings{ind};
      new_enc = [];
      if new_ind <= length(sys.state_encodings)
        new_enc = sys.state_encodings{new_ind};
      end
      sys.state_encodings{ind} = new_enc;
      sys.state_encodings{new_ind} = old_enc;
      mexBDD('set_state_enc', sys.BDD_system_ID, uint32(ind), uint32(new_enc));
      mexBDD('set_state_enc', sys.BDD_system_ID, uint32(new_ind), uint32(old_enc));
    end
    
    function add_progress_group(sys, U, G)
      mexBDD('add_pg', sys.BDD_system_ID, U, G);
    end
    
    function rm_progress_group(sys, ind)
      mexBDD('rm_pg', sys.BDD_system_ID, ind);
    end
    
    function [U, G] = get_progress_group(sys, ind)
      [U, G] = mexBDD('read_pg', sys.BDD_system_ID, ind);
      U = sort(U);
      G = sort(G);
    end
    
    function val = has_superior_pg(sys, U, G)
      val = mexBDD('check_sup_pg', sys.BDD_system_ID, U, G);
    end
    
    function inds = get_membership_pg(sys, G)
      inds = mexBDD('check_mem_in_G_pg', sys.BDD_system_ID, uint32(G));
    end
    
    function add_to_pg(sys, pg_inds, G)
      mexBDD('add_to_pg', sys.BDD_system_ID, pg_inds, G);
    end
    
    function [V, cont] = pre(sys, X, U, quant1, quant2)
      if nargout == 1
        V = mexBDD('pre', sys.BDD_system_ID, X, U, quant1, quant2);
      elseif nargout == 2
        [V, contID] = mexBDD('pre', sys.BDD_system_ID, X, U, quant1, quant2);
        cont = Controller.fromBDD(contID);
      end
      V = sort(V);
    end
    
    function [W, Cw, cont] = pre_pg(sys, V, B, quant1)
      if nargout == 3
        [W, Cw, contID] = mexBDD('pre_pg', sys.BDD_system_ID, V, B, quant1);
        Cw = sort(Cw);
        cont = Controller.fromBDD(contID);
      elseif nargout == 2
        [W, Cw] = mexBDD('pre_pg', sys.BDD_system_ID, V, B, quant1);
        Cw = sort(Cw);
      elseif nargout == 1
        W = mexBDD('pre_pg', sys.BDD_system_ID, V, B, quant1);
      end
      W = sort(W);
    end
    
    function [W, Cw, cont] = pginv(sys, U, G, Z, B, quant1)
      if nargout == 1
        W = mexBDD('pg_inv', sys.BDD_system_ID, U, G, Z, B, quant1);
      elseif nargout == 2
        [W, Cw] = mexBDD('pg_inv', sys.BDD_system_ID, U, G, Z, B, quant1);
        Cw = sort(Cw);
      elseif nargout == 3
        [W, Cw, contID] = mexBDD('pg_inv', sys.BDD_system_ID, U, G, Z, B, quant1);
        Cw = sort(Cw);
        cont = Controller.fromBDD(contID);
      end
      W = sort(W);
    end
    
    function [V, Cv, cont] = win_until(sys, B, P, quant1)
      if nargout == 1
        V = mexBDD('win_until', sys.BDD_system_ID, B, P, quant1);
      elseif nargout == 2
        [V, Cv] = mexBDD('win_until', sys.BDD_system_ID, B, P, quant1);
        Cv = sort(Cv);
      elseif nargout == 3
        [V, Cv, contID] = mexBDD('win_until', sys.BDD_system_ID, B, P, quant1);
        Cv = sort(Cv);
        cont = Controller.fromBDD(contID);
      end
      V = sort(V);
    end
    
    function [V, Cv, cont] = win_until_and_always(sys, A, B, P, quant1)
      if nargout == 1
        V = mexBDD('win_until_and_always', sys.BDD_system_ID, A, B, P, quant1);
      elseif nargout == 2
        [V, Cv] = mexBDD('win_until_and_always', sys.BDD_system_ID, A, B, P, quant1);
        Cv = sort(Cv);
      elseif nargout == 3
        [V, Cv, contID] = mexBDD('win_until_and_always', sys.BDD_system_ID, A, B, P, quant1);
        Cv = sort(Cv);
        cont = Controller.fromBDD(contID);
      end
      V = sort(V);
    end
    
    function [V, Cv, cont] = win_intermediate(sys, A, B, P, C_list, quant1)
      if nargout == 1
        V = mexBDD('win_intermediate', sys.BDD_system_ID, A, B, P, C_list, quant1);
      elseif nargout == 2
        [V, Cv] = mexBDD('win_intermediate', sys.BDD_system_ID, A, B, P, C_list, quant1);
        Cv = sort(Cv);
      elseif nargout == 3
        [V, Cv, contID] = mexBDD('win_intermediate', sys.BDD_system_ID, A, B, P, C_list, quant1);
        Cv = sort(Cv);
        cont = Controller.fromBDD(contID);
      end
      V = sort(V);
    end
    
    function [V, Cv, cont] = win_primal(sys, A, B, C_list, quant1, quant2, V)
      if nargin < 7
        V = [];
      end
      if nargout == 1
        V = mexBDD('win_primal', sys.BDD_system_ID, A, B, C_list, quant1, quant2, V);
      elseif nargout == 2
        [V, Cv] = mexBDD('win_primal', sys.BDD_system_ID, A, B, C_list, quant1, quant2, V);
        Cv = sort(Cv);
      elseif nargout == 3
        [V, Cv, contID] = mexBDD('win_primal', sys.BDD_system_ID, A, B, C_list, quant1, quant2, V);
        Cv = sort(Cv);
        cont = Controller.fromBDD(contID);
      end
      V = sort(V);
    end
    
    function time = win_primal_time(sys, A, B, C_list, quant1, quant2, V)
      if nargin < 7
        V = [];
      end
      time = mexBDD('win_primal_time', sys.BDD_system_ID, A, B, C_list, quant1, quant2, V);
    end
    
    function print_states(sys)
      mexBDD('print_all_states', sys.BDD_system_ID);
    end
    
    function print_transitions(sys)
      mexBDD('print_trans', sys.BDD_system_ID);
    end
    
    function ret = num_trans(sys)
      ret = mexBDD('num_trans', sys.BDD_system_ID);
    end
    
    function print_encs(sys)
      mexBDD('print_encs', sys.BDD_system_ID);
    end
    
    function debug(sys)
      mexBDD('debug', sys.BDD_system_ID);
    end
    
    function count = count_nodes(sys)
      count = mexBDD('count_nodes', sys.BDD_system_ID);
    end
    
    function count = count_dead_nodes(sys)
      count = mexBDD('count_dead_nodes', sys.BDD_system_ID);
    end
    
    function reorder(sys, bound)
      mexBDD('reorder', sys.BDD_system_ID, bound);
    end
    
    function dyn_reordering(sys, is_on)
      mexBDD('toggle_reorder', sys.BDD_system_ID, is_on);
    end
    
    function count = count_system_nodes(sys)
      count = mexBDD('count_system_nodes', sys.BDD_system_ID);
    end
    
    function read_var_order(sys)
      mexBDD('read_var_order', sys.BDD_system_ID);
    end
    
    function sobj = saveobj(obj)
      sobj = obj;
      % returns data fields to save (as cell array) and creates
      % temporary file temp_sys.dump
      save_data = mexBDD('save', obj.BDD_system_ID);
      
      % save data fields
      sobj.s_in_inds = save_data{1};
      sobj.s_out_inds = save_data{2};
      sobj.a_inds = save_data{3};
      sobj.pg_num = save_data{4};
      
      % read and save BDDs of system as string
      s = dir('temp_sys.dump');
      sobj.BDD_file_string = blanks(s.bytes);
      fileID = fopen('temp_sys.dump');
      pos = 1;
      while ~feof(fileID)
        line = fgets(fileID);
        line_len = length(line);
        sobj.BDD_file_string(pos:(pos+line_len-1)) = line;
        pos = pos + line_len;
      end
      fclose(fileID);
      delete temp_sys.dump
      if pos < length(obj.BDD_file_string)
        sobj.BDD_file_string = sobj.BDD_file_string(1:(pos-1));
      end
    end
      
  end
  
  methods (Static)
    function lobj = loadobj(obj)
      disp('hello');
      s.s_var_num = obj.s_var_num;
      s.a_var_num = obj.a_var_num;
      
      s.state_encodings = obj.state_encodings;
      s.action_encodings = obj.action_encodings;
      
      s.s_in_inds = obj.s_in_inds;
      obj.s_in_inds = [];
      s.s_out_inds = obj.s_out_inds;
      obj.s_out_inds = [];
      s.a_inds = obj.a_inds;
      obj.a_inds = [];
      s.pg_num = obj.pg_num;
      
      % create dump file for reading by dddmp
      dumpID = fopen('temp_sys.dump', 'w');
      fprintf(dumpID, '%s', obj.BDD_file_string);
      fclose(dumpID);
      obj.BDD_file_string = [];
      
      lobj = obj; 
      lobj.BDD_system_ID = mexBDD('load', s);
      delete temp_sys.dump
    end
  end
end

