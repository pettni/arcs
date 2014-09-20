classdef Partition<handle

	properties (SetAccess=protected)
		dim;		% Dimension of partition
		cell_list;	% List of cells
		domain;		% Domain
		adjacency;	% Sparse matrix containing adjacency information
		adjacency_outside;	% Which cells are at the boundary of the domain?
	end

	methods
		function r = Partition(domain)
			% PARTITION: Create a Partition object
			%
			% SYNTAX
			% ------
			%
			%	part = Partition(domain)
			% 
			% INPUT
			% -----
			%	
			%	domain 	domain of the Partition
			%
			% SEE ALSO
			% --------
			%
			% 	Partition/add_area (used to refine the partition), 
			% 	split_cell, get_neighbors, Partition/plot, get_all_aps, get_cells_with_ap, add_aps
			r.dim = domain.dim;
			r.domain = domain;
			r.cell_list = [domain];
			r.adjacency = [];
		end

		function add_area(part, area)
			% ADD_AREA: Adds a cell with to the partition. This new cell is subtracted from existing cells,
			% so existing AP's will not be preserved.
			%
			% SYNTAX
			% ------
			%
			%	part.add_area(cell)
			% 
			% INPUT
			% -----
			%	
			%	cell 	the new cell to add
			%
			% OUTPUT
			% ------
			%
			% 	none
			part.cell_list = [mldivide(part.cell_list, area) area];
		end

		function compute_adjacency(part)
			% Generate and store adjacency matrix
		    N = length(part.cell_list);
		    adjacency = sparse(N, N);
		    adjacency_outside = zeros(1,N);
		    for i=1:N
		        for j=i+1:N
		            [isn, d] = isNeighbor(part.cell_list(i), part.cell_list(j));
		            if isn
		                adjacency(i,j) = d;
		                adjacency(j,i) = d;
		            end
		        end
		       	isi = contains_strictly(part.cell_list(i), part.domain);
		        if ~isi
		            adjacency_outside(i) = 1;
        		end
        	end
    		part.adjacency = adjacency;
    		part.adjacency_outside = adjacency_outside;
        end
        
        function [adj, dim] = get_neighbors(p,ind)
        	% Return indices of neighbors of cell number ind
        	if isempty(p.adjacency)
        		p.compute_adjacency();
        	end
        	adj = find(p.adjacency(ind,:));
        	dim = p.adjacency(ind,adj);
        end

        function ret = length(obj)
        	% Overload length
        	ret = length(obj.cell_list);
        end

        function varargout = subsref(obj,S)
        	% Overload indexing operator to access cells directly
			switch S(1).type
			  	case '.'
			  		% call builtin method to access class properties
			    	[varargout{1:nargout}] = builtin('subsref', obj, S);
			  	case '()'
			  		% delegate to cell_list indexing
			  		varargout{1} = obj.cell_list(S.subs{:});
		   	end 
		end

		function check(part)
			% Check that the cells cover the domain, and vice versa
			if part.domain ~= part.cell_list
    			error('Invalid partition');
			end
		end

		function ret = get_all_aps(part)
			% Return a list of all APs present in the partition
			ret = [];
			for c = part.cell_list
				ret = union(ret, c.ap);
			end
		end

		function ret = get_cells_with_ap(part, ap)
			% Return indices of cells with a given set of APs
			% if ap = [], return cells without an AP
			ret = [];
			if isempty(ap)
				for i = 1:length(part.cell_list)
					if isempty(part.cell_list(i).ap);
						ret(end+1) = i;
					end
				end
			else
				for i = 1:length(part.cell_list)
					if all(ismember(ap, part.cell_list(i).ap))
						ret(end+1) = i;
					end
				end
			end
		end

		function add_aps(part, list, ap)
			% Add a given AP to a list of cells
			for rec=part.cell_list(list)
				rec.add_ap(ap);
			end
		end

    end
end