function [V, C, cont] = win_until_or_always(ts, B_list, Z, quant1)
  % win_until_or_always: compute winning set of      
  %     or_i ( B_i U Z || []B_i )
  % under (quant1, 'forall')-controllability

  num_subset = 2^length(B_list);

  V_list = cell(1, num_subset-1);
  for j_ind = 1:num_subset-1
    V_list{j_ind} = uint32(1:ts.n_s);
  end

  if nargout > 2
    K_list = cell(1, num_subset-1);
  end

  while true
    changed = false;
    for j_ind = num_subset-1:-1:1
      jb = dec2bin(j_ind, length(B_list)) - '0';
      B_cap = uint32(1:ts.n_s);
      for j=find(jb)
        % For all indices
        B_cap = intersect(B_cap, B_list{j});
      end

      % U_{k subset J} V_list{k}
      X_union = [];
      for k_ind = 1:j_ind
        kb = dec2bin(k_ind, length(B_list)) - '0';
        if all(or(not(kb), jb))  % is kb subset of jb?
          % For all subsets
          X_union = union(X_union, V_list{k_ind});
        end
      end

      if nargout > 2
        [preX, K] = ts.pre(X_union, [], quant1, 0);
        K.restrict_to(B_cap);
        K_list{num_subset-j_ind} = K;
      else
        preX = ts.pre(X_union, [], quant1, 0);
      end

      V_j = union(Z, intersect(B_cap, preX));
      if length(V_j) ~= length(V_list{j_ind})
        changed = true;
        V_list{j_ind} = V_j;
      end
    end
    if ~changed
      break
    end
  end

  V_log = false(1, ts.n_s);
  for j_ind = 1:num_subset-1
    V_log(V_list{j_ind}) = true;
  end
  V = zeros(1, sum(V_log), 'uint32');
  V(:) = find(V_log);

  C = [];  % todo

  if nargout > 2
    cont = Controller(V_list, K_list, 'reach', 'win_until_or_always');
  end
end
