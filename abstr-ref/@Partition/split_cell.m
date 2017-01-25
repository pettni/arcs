function [ind1, ind2] = split_cell(part, ind, dim)
  % SPLIT_CELL: Split a cell in the partition, which increases the number of cells by 1.
  % The adjacency matrix is updated automatically.
  % 
  % SYNTAX
  % ------
  %
  % [io,jo] = part.split_cell()
  % [io,jo] = part.split_cell(i)
  % [io,jo] = part.split_cell(i, dim)
  % 
  % INPUT
  % -----
  % 
  % i   number of cell to split (default: largest cell)
  %   dim   dimension to split along (default: largest dimension)
  %
  % OUTPUT
  % ------
  %
  %   [io, jo]  indices of new cells
  if nargin<2
    [~, ind] = max(volume(part.cell_list));
  end

  if nargin<3
    [~, dim] = max(part.cell_list(ind).xmax - part.cell_list(ind).xmin);
  end

  old_adj_ind = part.get_neighbors(ind);
  N = length(part);

  % split cell number ind along dimension dim
  [p1 p2] = split(part.cell_list(ind),dim);

  %%% contstruct adjacency to be filled %%%
  new_adj = [part.adjacency zeros(N,1); zeros(1, N+1)];
  new_adj(ind,:) = 0;
  new_adj(:,ind) = 0;

  % mutual adjacency along dim
  new_adj(ind, N+1) = dim;
  new_adj(N+1, ind) = dim;

  % determine kept adj
  for num_oldadj=old_adj_ind
    [isn, d] = isNeighbor(p1, part.cell_list(num_oldadj));
    if isn
      new_adj(ind, num_oldadj) = d;
      new_adj(num_oldadj, ind) = d;
    end
    [isn, d] = isNeighbor(p2, part.cell_list(num_oldadj));
    if isn
      new_adj(N+1, num_oldadj) = d;
      new_adj(num_oldadj, N+1) = d;
    end
  end

  new_adj_out = [part.adjacency_outside 0];
  new_adj_out(ind) = 0;
  % adjacent to outside
  if part.adjacency_outside(ind)
    isi = contains_strictly(p1, part.domain);
    if ~isi
      new_adj_out(ind) = 1;
    end
    isi = contains_strictly(p2, part.domain);
    if ~isi
      new_adj_out(N+1) = 1;
    end
  end

  % Add new cells and adjacency matrix
  part.cell_list(ind) = p1; % add one at old position
  part.cell_list(end+1) = p2;     % and one at the end
  part.adjacency = new_adj;
  part.adjacency_outside = new_adj_out;
  ind1=ind; ind2=N+1;
end
