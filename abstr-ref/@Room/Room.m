classdef Room<handle
    %ROOM Represents a room in a building structure
    %   An object used for object referensing and building construction
    
    properties
        dimensions          % [m] dimensions of room [x, y, z]
        A_fac = 0;          % [m^2] area to building exterior
        rho_a = 1.2041;     % [kg/m^3] density, air
        C                   % [J/K] heat capacitance, air
        c_a = 1012;         % [J/kgK] specific heat, air	
        r_o = 2.356*18/36;  % [m^2K/W] Room/outside r value
        r_c = 0.125         % [m^2K/W] Room/slab r value
        r_room = 0.79       % [m^2K/W] Room/room r value
        T_o = 30;           % [C] Outside temperature
        heat                % [J/m^2] added heat 
        R_vals = [Inf];     % [K/W] R value to other structures (slab first)
        neighs = [];        % neighbouring rooms
        neigh_a = [];       % common wall areas of neighbours
        slab                % present slab
        id                  % identity number of room
        pos                 % position of room in building
    end
    
    methods (Static)
        function r = get_room(type, pos, rotate)
            
            r = Room();
            if type == 1
                r.dimensions = [8,6,3];
                r.heat = 6;
                r.A_fac = 24;
            elseif type == 2
                r.dimensions = [6,6,3];
                r.heat = 8;
                r.A_fac = 18;
            else
                error('Unknown room type');
            end
            
            if nargin < 3
                rotate = false;
            end
            
            if rotate
                r.dimensions([2,1]) = r.dimensions([1,2]);
            end
            
            r.pos = pos;
%             for i = 1:3
%                 r.A_fac = r.A_fac + 2 * prod(r.dimensions(setdiff(1:3, i)));
%             end
%             r.A_fac = r.A_fac - r.floor_area();
            r.C = 5 * r.c_a * r.rho_a * prod(r.dimensions);
            Room.getset_room_num(1);
            r.id = Room.getset_room_num();
        end
        
        function num = getset_room_num(val_mod, clear)
            persistent room_num
            if isempty(room_num)
                room_num = 0;
            end
            
            if nargin > 1
                room_num = 0;
            end
            if nargin < 1
                num = room_num;
            else
                room_num = room_num + val_mod;
            end
        end
    end
    
    methods
        function r = Room()
        end
        
        function setSlab(r, slab)
%             if isempty(r.slab)
%                 r.A_fac = r.A_fac - r.floor_area();
%             end
            r.slab = slab;
            r.slab.total_area = r.slab.total_area + r.floor_area();
            slab.add_room(r);
        end
        
        function a = floor_area(r)
            a = prod(r.dimensions([1,2]));
        end
        
        function add_neigh(r, neigh)
            dist = r.pos - neigh.pos;
            if sum(abs(dist)) ~= 1
                error('invalid neighbour');
            end
            
            common_a = 0;
            if dist(1) ~= 0
                common_a = min([r.dimensions(2), neigh.dimensions(2)]) ...
                           * min([r.dimensions(3), neigh.dimensions(3)]);
            else
                common_a = min([r.dimensions(1), neigh.dimensions(1)]) ...
                           * min([r.dimensions(3), neigh.dimensions(3)]);
            end
%            r.A_fac = r.A_fac - common_a;
            
            if isempty(r.neighs)
                r.neighs = [neigh];
            else
                r.neighs(end+1) = neigh;
            end
            r.neigh_a(end+1) = common_a;
        end
        
        function [A_dyn, E_dyn, K_dyn] = get_dyn(r)
            slab_num = Slab.getset_slab_num();
            room_num = Room.getset_room_num();
            var_num = slab_num + room_num;
            place = r.id+slab_num;
            A_dyn = zeros(var_num);
            E_dyn = zeros(var_num, room_num);
            K_dyn = zeros(var_num, 1);
            
            % neighbouring rooms
            for i = 1:length(r.neighs)
                neigh = r.neighs(i);
                neigh_var = neigh.id + slab_num;
                R = r.r_room/r.neigh_a(i);
                A_dyn(place, neigh_var) = 1/R;
                A_dyn(place, place) = A_dyn(place,place) - 1/R;
            end
            
            % slab
            R_c = r.r_c / r.floor_area();
            if ~isempty(r.slab)
                A_dyn(place, r.slab.id) = 1/R_c;
                A_dyn(place,place) = A_dyn(place,place) - 1/R_c;
            end
            
            % exterior
            R_o = r.r_o / r.A_fac;
            A_dyn(place,place) = A_dyn(place,place) - 1/R_o;
            K_dyn(place) = K_dyn(place) + r.T_o/R_o;
            
            % heat
            K_dyn(place) = K_dyn(place) + r.heat*r.floor_area();
            
            % disturbance
            E_dyn(place, r.id) = r.floor_area;
            
            A_dyn = sparse(A_dyn)/r.C;
            E_dyn = E_dyn/r.C;
            K_dyn = sparse(K_dyn/r.C);
            
        end
        
        
    end
end

