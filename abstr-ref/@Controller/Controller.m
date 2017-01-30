classdef Controller<handle

  properties (SetAccess=protected)
    sets;           % List of sets where subcontrollers are active
    subcontrollers; % List of Controller or containers.Map
    control_type;   % One of 'reach', 'recurrance', and 'simple'
    mem_var;        % Memory variable needed for higher-level controller
  end

  methods
    function cont = Controller(set_list, c_list, c_type)
      % Some checks
      if length(set_list) ~= length(c_list)
        error('set_list and c_list must be same size')
      end

      cont.sets = set_list;
      cont.subcontrollers = c_list;
      cont.control_type = c_type;
      cont.mem_var = [];
    end

    function varargout = subsref(obj,S)
      % Overload indexing operator to access cells directly
      switch S(1).type
        % call builtin method to access class properties
      case '.' 
        [varargout{1:nargout}] = builtin('subsref', obj, S);
        % delegate to cell_list indexing
      case '()'
        varargout{1} = obj.get_input(S.subs{:});
      end 
    end

    function  a_list = get_input(cont, state)
      % Return controller input for 'state'
      if cont.control_type == 'simple'
        a_list = cont.subcontrollers(state);
      end
    end

    function restrict_to(cont, r_set)
      if cont.control_type == 'simple'
        curr_keys = keys(cont.subcontrollers);
        cont.subcontrollers.remove(setdiff(cell2mat(curr_keys), r_set));
      end
    end
  end
end