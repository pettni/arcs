function [ind1, ind2] = split_cell(part, ind, dim)
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
  [p1, p2] = split(part.cell_list(ind),dim);

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

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % IF THERE IS A TS, UPDATE IT %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if isempty(part.ts)
    return
  end
  
  if strcmp(part.ts.sys_setting, TransSyst.sparse_set)
    % Store transitions pertaining to split state
    % to reduce computation
    trans_in = [];
    trans_in_act = [];
    trans_out = [];
    trans_out_act = [];
    outdomain_act = [];
    non_transient_act = [];

    % remove all transitions pertaining to ind1
    % and store information
%     c = length(part.ts.action);
    for i=part.ts.num_trans():-1:1
      if part.ts.state1(i) == ind1 && part.ts.state2(i) == part.ts.n_s
        outdomain_act(end+1) = part.ts.action(i);
      elseif part.ts.state1(i) == ind1 && part.ts.state2(i) == ind1
        non_transient_act(end+1) = part.ts.action(i);
      elseif part.ts.state1(i) == ind1 
        trans_out(end+1) = part.ts.state2(i);
        trans_out_act(end+1) = part.ts.action(i);
      elseif part.ts.state2(i) == ind1 
        trans_in(end+1) = part.ts.state1(i);
        trans_in_act(end+1) = part.ts.action(i);
      end
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

    % increase state counter
    part.ts.n_s = part.ts.n_s + 1;
    part.ts.array_computed = false;

    trans_to_add = zeros(0, 3);
    % Transitions into union
    for i = 1:length(trans_in)
      in_state = trans_in(i);
      in_act = trans_in_act(i);
      if is_trans(part.cell_list(in_state), ...
                  part.cell_list(ind1), ...
                  part.dyn_list{in_act}, part.d_rec)
        part.ts.add_transition(in_state, ind1, in_act);
      end
      if is_trans(part.cell_list(in_state), ...
                  part.cell_list(ind2), ...
                  part.dyn_list{in_act}, part.d_rec)
        part.ts.add_transition(in_state, ind2, in_act);
      end
    end

    % Transitions out of union
    for i = 1:length(trans_out)
      out_state = trans_out(i);
      out_act = trans_out_act(i);
      if is_trans(part.cell_list(ind1), ...
                  part.cell_list(out_state), ...
                  part.dyn_list{out_act}, part.d_rec)
        part.ts.add_transition(ind1, out_state, out_act);
      end
      if is_trans(part.cell_list(ind2), ...
                  part.cell_list(out_state), ...
                  part.dyn_list{out_act}, part.d_rec)
        part.ts.add_transition(ind2, out_state, out_act);
      end
    end

    % Transitions between the two cells
    for act_num = 1:length(part.dyn_list)
      if is_trans(part.cell_list(ind1), ...
                  part.cell_list(ind2), ...
                  part.dyn_list{act_num}, part.d_rec)
        part.ts.add_transition(ind1, ind2, act_num);
      end
      if is_trans(part.cell_list(ind2), ...
                  part.cell_list(ind1), ...
                  part.dyn_list{act_num}, part.d_rec)
        part.ts.add_transition(ind2, ind1, act_num);
      end
    end

    % Out-of-domain
    for i = 1:length(outdomain_act)
      out_act = outdomain_act(i);
      if is_trans_out(part.cell_list(ind1), ...
                      part.domain, ...
                      part.dyn_list{out_act}, part.d_rec)
        part.ts.add_transition(ind1, N+2, out_act);
      end
      if is_trans_out(part.cell_list(ind2), ...
                      part.domain, ...
                      part.dyn_list{out_act}, part.d_rec)
        part.ts.add_transition(ind2, N+2, out_act);
      end
    end

    % Self transitions
    for i = 1:length(non_transient_act)
      nt_act = non_transient_act(i);
      if ~is_transient(part.cell_list(ind1), ...
                       {part.dyn_list{nt_act}}, part.d_rec)
        part.ts.add_transition(ind1, ind1, nt_act);
      end
      if ~is_transient(part.cell_list(ind2), ...
                       {part.dyn_list{nt_act}}, part.d_rec)
        part.ts.add_transition(ind2, ind2, nt_act);
      end
    end
    % update existing progress groups
    for i = 1:length(part.ts.pg_G)
      if ismember(ind1, part.ts.pg_G{i})
        part.ts.pg_G{i} = union(part.ts.pg_G{i}, ind2);
      end
    end

    % Check for new progress groups in transient regions
    for i=1:length(part.trans_reg_U)
      if (part.trans_reg_rec{i}.contains(part.cell_list(ind1)) || ...
          part.trans_reg_rec{i}.contains(part.cell_list(ind2)))
        G_cand = part.get_cells_inside(part.trans_reg_rec{i});
        U_cand = part.trans_reg_U{i};
        if ~part.ts.has_superior_pg(U_cand, G_cand)
          part.ts.add_progress_group(U_cand, G_cand);
        end
      end
    end
  elseif strcmp(part.ts.sys_setting, TransSyst.bdd_set)
    % Store transitions pertaining to split state
    % to reduce computation
    trans_in = [];
    trans_in_act = [];
    trans_out = [];
    trans_out_act = [];
    outdomain_act = [];
    non_transient_act = [];

    % get all transitions pertaining to ind1
    [state1, action, state2] = part.ts.bdd_sys.get_trans_with_state(ind1);
    trans_num = length(state1); % number of transitions to investigate

    % remove all transitions pertaining to ind1
    % and store information
    for i=trans_num:-1:1
      if state1(i) == ind1 && state2(i) == part.ts.n_s
        outdomain_act(end+1) = action(i);
      elseif state1(i) == ind1 && state2(i) == ind1
        non_transient_act(end+1) = action(i);
      elseif state1(i) == ind1 
        trans_out(end+1) = state2(i);
        trans_out_act(end+1) = action(i);
      elseif state2(i) == ind1 
        trans_in(end+1) = state1(i);
        trans_in_act(end+1) = action(i);
      end
    end
    
    part.ts.bdd_sys.remove_trans_with_state(ind1);

    % move last (outside) state forward---should only have 
    % outgoing transitions
    % Move out state down in state list
    % increase state counter and update encodings
    part.ts.add_state(ind1);

    part.ts.bdd_sys.change_state(part.ts.n_s-1, part.ts.n_s);

    trans_to_add = zeros(0, 3);
    % Transitions into union
    for i = 1:length(trans_in)
      in_state = trans_in(i);
      in_act = trans_in_act(i);
      if is_trans(part.cell_list(in_state), ...
                  part.cell_list(ind1), ...
                  part.dyn_list{in_act}, part.d_rec)
        trans_to_add(end+1, :) = [in_state, ind1, in_act];
      end
      if is_trans(part.cell_list(in_state), ...
                  part.cell_list(ind2), ...
                  part.dyn_list{in_act}, part.d_rec)
        trans_to_add(end+1, :) = [in_state, ind2, in_act];
      end
    end

    % Transitions out of union
    for i = 1:length(trans_out)
      out_state = trans_out(i);
      out_act = trans_out_act(i);
      if is_trans(part.cell_list(ind1), ...
                  part.cell_list(out_state), ...
                  part.dyn_list{out_act}, part.d_rec)
        trans_to_add(end+1,:) = [ind1, out_state, out_act];
      end
      if is_trans(part.cell_list(ind2), ...
                  part.cell_list(out_state), ...
                  part.dyn_list{out_act}, part.d_rec)
        trans_to_add(end+1, :) = [ind2, out_state, out_act];
      end
    end

    % Transitions between the two cells
    for act_num = 1:length(part.dyn_list)
      if is_trans(part.cell_list(ind1), ...
                  part.cell_list(ind2), ...
                  part.dyn_list{act_num}, part.d_rec)
        trans_to_add(end+1, :) = [ind1, ind2, act_num];
      end
      if is_trans(part.cell_list(ind2), ...
                  part.cell_list(ind1), ...
                  part.dyn_list{act_num}, part.d_rec)
        trans_to_add(end+1, :) = [ind2, ind1, act_num];
      end
    end

    % Out-of-domain
    for i = 1:length(outdomain_act)
      out_act = outdomain_act(i);
      if is_trans_out(part.cell_list(ind1), ...
                      part.domain, ...
                      part.dyn_list{out_act}, part.d_rec)
        trans_to_add(end+1, :) = [ind1, N+2, out_act];
      end
      if is_trans_out(part.cell_list(ind2), ...
                      part.domain, ...
                      part.dyn_list{out_act}, part.d_rec)
        trans_to_add(end+1, :) = [ind2, N+2, out_act];
      end
    end

    % Self transitions
    for i = 1:length(non_transient_act)
      nt_act = non_transient_act(i);
      if ~is_transient(part.cell_list(ind1), ...
                       {part.dyn_list{nt_act}}, part.d_rec)
        trans_to_add(end+1, :) = [ind1, ind1, nt_act];
      end
      if ~is_transient(part.cell_list(ind2), ...
                       {part.dyn_list{nt_act}}, part.d_rec)
        trans_to_add(end+1, :) = [ind2, ind2, nt_act];
      end
    end

    part.ts.add_transition(trans_to_add(:, 1), trans_to_add(:, 2), trans_to_add(:, 3));

    % update existing progress groups

    pg_indices = part.ts.bdd_sys.get_membership_pg(ind1);
    part.ts.bdd_sys.add_to_pg(pg_indices, ind2);

    % Check for new progress groups in transient regions
    for i=1:length(part.trans_reg_U)
      if (part.trans_reg_rec{i}.contains(part.cell_list(ind1)) || ...
          part.trans_reg_rec{i}.contains(part.cell_list(ind2)))
        G_cand = part.get_cells_inside(part.trans_reg_rec{i});
        U_cand = part.trans_reg_U{i};
        if ~part.ts.has_superior_pg(U_cand, G_cand)
          part.ts.add_progress_group(U_cand, G_cand);
        end
      end
    end
end
