function [] = search_trans_reg(part, search_depth)
  % search_trans_reg: search for transient regions 
  % in partition domain

  search_set_queue = [part.domain];
  search_set_mindim = [1];
  search_set_numsplit = [0];

  while length(search_set_queue) > 0
    % Take from queue
    search_set = search_set_queue(1);
    min_dim = search_set_mindim(1);
    num_split = search_set_numsplit(1);

    split_cells = [];
    split_mindim = [];
    split_numsplit = [];

    if num_split+1 <= search_depth
      % Add smaller sets
      for j=min_dim:search_set.dim
        [s1 s2] = split(search_set, j);
        split_cells = [split_cells s1 s2];
        split_mindim = [split_mindim j j];
        split_numsplit = [split_numsplit num_split+1 num_split+1];
      end
    end

    % Update queue
    search_set_queue = [search_set_queue(2:end) split_cells];
    search_set_mindim = [search_set_mindim(2:end) split_mindim];
    search_set_numsplit = [search_set_numsplit(2:end), split_numsplit];

    U = 1:part.ts.n_a;
    for U_size = length(U):-1:1
      for Up = nchoosek(U, U_size)'
        if ~part.has_superior_trans_reg(Up, search_set)
          if is_transient(search_set, part.dyn_list(Up'), part.d_rec)
            part.trans_reg_U{end+1} = Up;
            part.trans_reg_rec{end+1} = search_set;
          end 
        end
      end
    end    
  end
end