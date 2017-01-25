classdef Partition<handle

  properties (SetAccess=protected)
    dim;  % Dimension of partition
    cell_list; % List of cells
    domain;  % Domain
    adjacency; % Sparse matrix containing adjacency information
    adjacency_outside; % Which cells are at the boundary of the domain?
  end

  methods
    function r = Partition(domain)
      % PARTITION: Create a Partition object
      r.dim = domain.dim;
      r.domain = domain;
      r.cell_list = [domain];
      r.adjacency = [];
    end

    function add_area(part, area)
      % Adds a cell with to the partition. This new cell
      % is subtracted from existing cells, so existing AP's
      % may not be preserved.
      part.cell_list = [mldivide(part.cell_list, area) area];
    end

    function ret = find_cell(part, point)
      % find cell that contains given point
      ret = -1;
      for i=1:length(part.cell_list)
        if isInside(part.cell_list(i), point)
          ret = i;
          return
        end
      end
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
        % call builtin method to access class properties
      case '.' [varargout{1:nargout}] = builtin('subsref', obj, S);
        % delegate to cell_list indexing
      case '()' varargout{1} = obj.cell_list(S.subs{:});
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
      ret ={};
      for c = part.cell_list
        ret = union(ret, c.ap);
      end
    end

    function ret = get_cells_with_ap(part, ap)
      % Return indices of cells with a given set of APs
      % if ap = {}, return cells without an AP
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