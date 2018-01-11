classdef Slab<handle
    %Slab Represents a concrete slab added to floor of rooms in which
    %   chilling water flows.
    %   Simply an object which represents a slab. Is used to provide object
    %   object referencing when constructing building and dynamics.
    
    properties
        thickness = 0.25;      % [m] thickness of slab
        rho = 2300;         % [kg/m3] density, concrete
        c = 750;			% [J/kgK] specific heat, concrete
        r_w = 0.102;        % [m^2K/W] Water/slab r value
        T_w = 18;           % [C] water temp.
        total_area = 0;     % [m^2] total area of slab
        rooms = {};         % neighbouring rooms
        id;
    end
    
    methods(Static)
        function s = get_slab()
            s = Slab();
            Slab.getset_slab_num(1);
            s.id = Slab.getset_slab_num();
        end
        
        function num = getset_slab_num(val, clear)
            persistent slab_num
            if isempty(slab_num)
                slab_num = 0;
            end
            if nargin > 1
                slab_num = 0;
            end
            if nargin < 1
                num = slab_num;
            else
                slab_num = slab_num + val;
            end
        end
    end
    
    methods
        function s = Slab()
        end
        
        function v = volume(s)
            % calculate volume of slab
            v = s.thickness * s.total_area;
        end
        
        function C = heat_cap(s)
            C = s.c * s.rho * volume(s);
        end
        
        function add_room(s, room)
            if isempty(s.rooms)
                s.rooms = [room];
            elseif ~ismember(room.id, [s.rooms.id])
                s.rooms(end+1) = room;
            else
                error(['slab ', num2str(s.id), ' already in room ', num2str(room.id)]);
            end
        end
        
        function [A_dyn, K_dyn] = get_dyn(s, is_on)
            slab_num = Slab.getset_slab_num();
            room_num = Room.getset_room_num();
            var_num = slab_num + room_num;
            place = s.id;
            A_dyn = zeros(var_num);
            K_dyn = zeros(var_num, 1);
            for i = 1:length(s.rooms)
                room = s.rooms(i);
                room_var = room.id + slab_num;
                R = room.r_c / room.floor_area();
                A_dyn(place,room_var) = 1/R;
                A_dyn(place,place) = A_dyn(place,place) - 1/R;
            end
            if is_on
                R_w = s.r_w/s.total_area;
                A_dyn(place, place) = A_dyn(place,place) - 1/R_w;
                K_dyn = sparse(zeros(var_num,1));
                K_dyn(place) = s.T_w/R_w;
            end
            
            A_dyn = sparse(A_dyn)/s.heat_cap();
            K_dyn = sparse(K_dyn)/s.heat_cap();
        end
        
        
    end
    
    
end

