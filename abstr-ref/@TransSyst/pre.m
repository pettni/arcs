function [V, cont] = pre(ts, X, U, quant1, quant2)
  % Compute pre(X) under (quant1, quant2)-controllability
  % and action set U. If U == [] all actions are considered
  %
  % Note: X must be sorted!
  % Returns a sorted set
  %
  % quant = false: forall  (non controllable)
  % quant = true: exists   (controllable)

  if strcmp(ts.sys_setting, TransSyst.bdd_set)
    if (nargout == 1)
      V = ts.bdd_sys.pre(X, U, quant1, quant2);
    elseif nargout == 2
      [V, cont] = ts.bdd_sys.pre(X,U,quant1,quant2);
    end
    return;
  end
  
  if nargout > 1
    Kmap = containers.Map('KeyType', 'uint32', 'ValueType', 'any');
  end

  if isempty(U)
    U = uint32(1:ts.n_a);
  end

  log_idx = false(1, ts.n_s);
  if ts.fast_enabled
    for i=1:length(X)
      for j = ts.fast_pre_all{X(i)}
        log_idx(j) = 1;
      end
    end
  else
    for i=1:ts.num_trans()
      if ismember(ts.state2(i), X) && ismember(ts.action(i), U)
        log_idx(ts.state1(i)) = 1;
      end
    end
  end
  
  for q = 1:ts.n_s
    if ~log_idx(q)
        continue
    end
    act_list = false(1, length(U));   % outcome per action 
    for i = 1:length(U)
      a = U(i);
      if ts.fast_enabled
        aPost = ts.fast_post{(a-1) * ts.n_s + q};
      else
        aPost = ts.post(q, a);
      end
      if quant2
        act_list(i) = any(builtin('_ismemberhelper',aPost, X));
      else
        act_list(i) = all(builtin('_ismemberhelper',aPost, X)) && ~isempty(aPost);
      end
    end
    if quant1
      if ~any(act_list)
        log_idx(q) = 0;
      elseif nargout > 1
        Kmap(q) = U(act_list); 
      end
    else
      if ~all(act_list)
        log_idx(q) = 0;
      elseif nargout > 1
        Kmap(q) = U;
      end
    end
  end

  V = zeros(1, sum(log_idx), 'uint32');
  V(:) = find(log_idx);
  if nargout > 1
    cont = Controller(V, Kmap, 'simple', 'pre');
  end
end
