function split_cell(part, ind, dim)
  % SPLIT_CELL: Split a cell in the partition, which increases the number of cells by 1.
  % The adjacency matrix is updated automatically.
  %
  % If there is an associated transition system, this is also updated
  % 
  % SYNTAX
  % ------
  %
  % part.split_cell()
  % part.split_cell(i)
  % part.split_cell(i, dim)
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

  %%% construct adjacency to be filled %%%
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

  if isempty(part.ts)
    return
  end

  % remove all transitions pertaining to ind1
  for i=part.ts.num_trans():-1:1
    if part.ts.state1(i) == ind1 || part.ts.state2(i) == ind1
      part.ts.state1(i) = [];
      part.ts.state2(i) = [];
      part.ts.action(i) = [];
    end
  end

  % move last (outside) state forward---should only have 
  % outgoing transitions
  for i=1:part.ts.num_trans()
    if part.ts.state2(i) == part.ts.n_s
      part.ts.state2(i) = part.ts.n_s+1;
    end
  end

  % update progress groups
  for i = 1:length(part.ts.pg_G)
    if ismember(ind1, part.ts.pg_G{i})
      part.ts.pg_G{i} = union(part.ts.pg_G{i}, part.ts.n_s);
    end
  end

  % increase state counter
  part.ts.n_s = part.ts.n_s + 1;
  part.ts.fast_enabled = false;

  % Re-establish transitions
  for act_num = 1:length(part.act_list)

    act = part.act_list{act_num}{1};
    dyn_type = part.act_list{act_num}{2};
    disturbance = part.act_list{act_num}{3};

    if strcmp(dyn_type, 'linear') && ~disturbance
      trans_fun = @(p1, p2) isTransLin(p1, p2, act);
      trans_out_fun = @(p1) isTransOutLin(p1, part.domain, act);
      transient_fun = @(p1) isTransientLin(p1, act);
    end

    for i = [ind1, ind2]
      for j = part.get_neighbors(i)
        if trans_fun(part.cell_list(i), part.cell_list(j))
          part.ts.add_transition(i, j, act_num);
        end
        if trans_fun(part.cell_list(j), part.cell_list(i))
          part.ts.add_transition(j, i, act_num);
        end
      end

      % Out-of-domain
      if trans_out_fun(part.cell_list(i))
        part.ts.add_transition(i, N+2, act_num);
      end

      % Self transitions
      if ~transient_fun(part.cell_list(i))
        part.ts.add_transition(i, i, act_num);
      end
    end
  end

  % TODO: check for new progress groups

end
