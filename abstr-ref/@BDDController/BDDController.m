classdef BDDController < handle
  %BDDCONTROLLER Summary of this class goes here
  %   Detailed explanation goes here
  
  properties (Access=private)
    BDD_cont_ID = -1;
  end
  
  methods
    function cont = BDDController(ID)
      cont.BDD_cont_ID = ID;
    end
    
    function restrict_to(cont, set)
      mexBDD('cont_restrict', cont.BDD_cont_ID, set);
    end
    
    function inputs = get_input(cont, states)
      inputs = mexBDD('cont_get_input', cont.BDD_cont_ID, states);
      inputs = sort(inputs);
    end
    
    function delete(cont)
      mexBDD('cont_delete', cont.BDD_cont_ID);
    end
  end
  
end

