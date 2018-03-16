classdef Building<handle
    %BUILDING Configuration of rooms and slabs
    
    properties
        rooms;
        slabs;
    end
    
    methods(Static)
        function build = create_building(room_pos, room_slabs, room_types)
          % create_building  Create object representing building, capable
          %                  creating building dynamics
          %
          %   building = create_building(room_pos, room_slabs, room_types)
          %   Creates a building with the rooms at 2D positions given by
          %   the rows of room_pos, where room_slabs(i) = k means room i 
          %   contains cooling slab k (k = 0 => no slab). The type of each
          %   room (type of dimension and other params) is specified by
          %   room_types, type being either 1 or 2.
          %
          build = Building(size(room_pos, 1), max(room_slabs));
            slabs = {};
            for i = 1:max(room_slabs)
                slabs{end+1} = Slab.get_slab();
            end
            for i = 1:size(room_pos,1)
                id = build.add_room(room_types(i), room_pos(i,:));
                if room_slabs(i) > 0
                    build.add_slab_to(id, slabs{room_slabs(i)});
                end
            end
            build.figure_out_neighs();
        end
    end
    
    methods
        function b = Building(room_num, slab_num)
            rooms(room_num) = Room;
            if slab_num > 0
                slabs(slab_num) = Slab;
                b.slabs = slabs;
            end
            b.rooms = rooms;
        end
        
        function delete(b)
            Room.getset_room_num(0, true);
            Slab.getset_slab_num(0, true);
        end
        
        function r_id = add_room(b, type, pos)
            room = Room.get_room(type, pos);
            r_id = room.id;
            b.rooms(r_id) = room;
        end
        
        function s_id = add_slab_to(b, room_id, slab)
            b.rooms(room_id).setSlab(slab);
            b.slabs(slab.id) = slab;
        end
        
        function figure_out_neighs(b)
            for i = 1:length(b.rooms)
                for j = (i+1):length(b.rooms)
                    if sum(abs(b.rooms(i).pos - b.rooms(j).pos)) == 1
                        b.rooms(i).add_neigh(b.rooms(j));
                        b.rooms(j).add_neigh(b.rooms(i));
                    end
                end
            end
        end
        
        function [A_dyn, E_dyn, K_dyn] = get_dyn(b)
            slab_num = Slab.getset_slab_num();
            room_num = Room.getset_room_num();
            var_num = slab_num + room_num;
            A_dyn = {};
            E_dyn = {};
            K_dyn = {};
            A_room = {};
            E_room = {};
            K_room = {};
            for i = 1:length(b.rooms)
                [A_room{i}, E_room{i}, K_room{i}] = b.rooms(i).get_dyn();
            end
            for m = (2^slab_num-1):-1:0
                config = de2bi(m);
                config = fliplr([config, zeros(1,slab_num - length(config))]);
                A = zeros(var_num);
                E = zeros(var_num, room_num);
                K = zeros(var_num, 1);
                
                for i = 1:length(b.rooms)
                    A = A + A_room{i};
                    E = E + E_room{i};
                    K = K + K_room{i};
                end
                
                for i = 1:length(b.slabs)
                    [A_slab, K_slab] = b.slabs(i).get_dyn(config(i));
                    A = A + A_slab;
                    K = K + K_slab;
                end
                
                A_dyn{m+1} = A;
                E_dyn{m+1} = E;
                K_dyn{m+1} = K;
            end
            
            A_dyn = fliplr(A_dyn);
            E_dyn = fliplr(E_dyn);
            K_dyn = fliplr(K_dyn);
        end
    end
    
end

