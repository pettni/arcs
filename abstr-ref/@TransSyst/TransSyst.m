classdef TransSyst<handle

  properties (SetAccess={?Partition})
    n_s;
    n_a;
    
    % BDD settings (sys = {'sparse', 'bdd'})
    sys_setting = '';
    
    % BDD system
    bdd_sys = [];
    
    % Transitions (sparse)
    state1;
    state2; 
    action;

    % Cell with one sparse transition matrix for each action
    % compute with trans_array_enable() when TS is fully constructed
    array_computed = false;
    trans_array = [];
    
    % Progress groups
    pg_U = {};
    pg_G = {};

    % Pre-computed pre_all and post maps for speed (sparse)
  end
  
  properties (Constant)
    % Setting strings
    sparse_set = 'sparse';
    bdd_set = 'bdd';
  end

  properties
    % Disable progress groups
    b_disable_pg = false;

    % Debugging flag
    b_debug = false;
  end

  methods
    function ts = TransSyst(n_s, n_a, setting, encoding_setting)
      ts.n_s = uint32(n_s);
      ts.n_a = uint32(n_a);
      if nargin < 3
        ts.sys_setting = TransSyst.sparse_set;
      else
        ts.sys_setting = setting;
      end
      if nargin < 3 || strcmp(setting, TransSyst.sparse_set)
        % Create a sparse TransSyst with n_s states and n_a actions
        ts.state1 = uint32([]);
        ts.state2 = uint32([]);
        ts.action = uint32([]);
        ts.trans_array = [];
      elseif strcmp(setting, TransSyst.bdd_set)
        ts.bdd_sys = BDDSystem(n_s, n_a, encoding_setting);
      else
        disp('Invalid TransSyst initialization');
      end
    end

    function numact = add_action(ts)
      % Add an action, return its number
      ts.n_a = ts.n_a + 1;
      numact = ts.n_a;
      ts.array_computed = false;
      
      if strcmp(ts.sys_setting, ts.bdd_set)
        ts.bdd_sys.add_action(ts.n_a);
      end
    end
    
    function numstate = add_state(ts, old_state)
      % Add a state, return its number
      ts.n_s = ts.n_s + 1;
      numstate = ts.n_s;
      ts.array_computed = false;
      
      if nargin > 1 && strcmp(ts.sys_setting, ts.bdd_set)
        ts.bdd_sys.add_state(old_state, ts.n_s);
      end
    end

    function add_transition(ts, s1, s2, a)
      % Add transition from s1 to s2 under action a

      if ts.b_debug
        assert(1 <= s1 && s1 <= ts.n_s)
        assert(1 <= s2 && s2 <= ts.n_s)
        assert(1 <= a && a <= ts.n_a)
      end
      
      if strcmp(ts.sys_setting, TransSyst.sparse_set)
        ts.state1 = [ts.state1; s1];
        ts.state2 = [ts.state2; s2];
        ts.action = [ts.action; a];
        ts.array_computed = false;
      elseif strcmp(ts.sys_setting, TransSyst.bdd_set)
        ts.bdd_sys.add_transition(s1, a, s2);
      end
    end

    function ret = has_superior_pg(ts, U, G)
      % Return true if a superior progress group
      % is already present
      
      if strcmp(ts.sys_setting, TransSyst.sparse_set)
        ret = false;
        for i=1:length(ts.pg_U)
          if all(ismember(U, ts.pg_U{i})) && ...
              all(ismember(G, ts.pg_G{i}))
            ret = true;
            return
          end
        end
      elseif strcmp(ts.sys_setting, TransSyst.bdd_set)
        ret = ts.bdd_sys.has_superior_pg(U, G);
      end
    end

    function add_progress_group(ts, U, G)
      % Add progress group G under modes U

      if ts.b_debug
        assert(all(1 <= U) && all(U <= ts.n_a))
        assert(all(1 <= G) && all(G <= ts.n_s))
      end

      if strcmp(ts.sys_setting, TransSyst.sparse_set)
        % Remove any pgs that are inferior
        for i=length(ts.pg_U):-1:1
          if all(ismember(ts.pg_U{i}, U)) && ...
              all(ismember(ts.pg_G{i}, G))
            ts.pg_U(i) = [];
            ts.pg_G(i) = [];
          end
        end

        ts.pg_U{end+1} = U;
        ts.pg_G{end+1} = G;
      elseif strcmp(ts.sys_setting, TransSyst.bdd_set)
        % Remove any pgs that are inferior
        pg_num = ts.bdd_sys.get_pg_num();
        for i=pg_num:-1:1
          [pg_U, pg_G] = ts.bdd_sys.get_progress_group(i);
          if all(ismember(pg_U, U)) && ...
              all(ismember(pg_G, G))
            ts.bdd_sys.rm_progress_group(i);
          end
        end
        
        ts.bdd_sys.add_progress_group(U, G);
      end
    end

    function ret = num_trans(ts)
      % Number of transitions
      if strcmp(ts.sys_setting, TransSyst.bdd_set)
        ret = ts.bdd_sys.num_trans();
      else
        ret = length(ts.state1);
      end
    end
    
    function trans_array_enable(ts)
      % Compute the array representation from list representation
      if ts.array_computed
        return
      end
      % build trans_array
      ts.trans_array = cell(ts.n_a,1);
      for i = 1:ts.n_a
          % change sub to idx s.t. fast to assign 1's
          idx_m = sub2ind([ts.n_s,ts.n_s], ts.state1(ts.action==i), ts.state2(ts.action==i));
          ts.trans_array{i} = sparse([], [], [], double(ts.n_s), double(ts.n_s), length(idx_m));
          ts.trans_array{i}(idx_m)=true;
      end
      ts.array_computed = true;
    end 
   
    %% BDD API
    function initializeBDD(ts)
      % Initialize the BDD in c
      
      % get initial number of variables
      for i = 1:ts.n_s
        ts.s_var_num = max([ts.s_var_num, length(ts.state_encodings{i})]);
      end
      for i = 1:ts.n_a
        ts.a_var_num = max([ts.a_var_num, length(ts.action_encodings{i})]);
      end
      
      manageBDD('initialize', ts.n_s, ts.s_var_num, ts.state_encodings, ts.a_var_num, ts.n_a, ts.action_encodings);
    end
    
    function toggle_setting(ts)
      % Converts the current Transition system to use the list
      % representation from a BDD representation
      % (given it has that setting)
      
      % get all transitions
      if strcmp(ts.sys_setting, TransSyst.bdd_set)
        [s_in, a, s_out] = ts.bdd_sys.get_trans_with_state(1:ts.n_s);
        ts.state1 = s_in;
        ts.action = a;
        ts.state2 = s_out;

        pg_num = ts.bdd_sys.get_pg_num();
        for i = 1:pg_num
          [ts.pg_U{i}, ts.pg_G{i}] = ts.bdd_sys.get_progress_group(i);
        end

        ts.array_computed = false;
        ts.sys_setting = TransSyst.sparse_set;
        ts.bdd_sys = [];
      end
    end
  end
end
