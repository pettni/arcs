function [ret, K] = pre(ts, X, U, quant1, quant2)
  % Compute pre(X) under (quant1, quant2)-controllability
  % and action set U. 
  %
  % Note: X must be sorted!
  % Returns a sorted set

  if nargout > 1
    K = containers.Map('KeyType', 'uint32', 'ValueType', 'any');
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
    act_list = zeros(1, ts.n_a);   % outcome per action 
    for a = 1:ts.n_a
        if ts.fast_enabled
          aPost = ts.fast_post{(a-1) * ts.n_s + q};
        else
          aPost = ts.post(q, a);
        end
        if strcmp(quant2, 'exists')
          act_list(a) = any(builtin('_ismemberhelper',aPost, X));
        else
          act_list(a) = all(builtin('_ismemberhelper',aPost, X)) && ~isempty(aPost);
        end
    end
    if strcmp(quant1, 'exists')
      if ~any(act_list)
        log_idx(q) = 0;
      elseif nargout > 1
        K(q) = find(act_list); 
      end
    else
      if ~all(act_list)
        log_idx(q) = 0;
      elseif nargout > 1
        K(q) = 1:ts.n_a;
      end
    end
  end

  ret = zeros(1, sum(log_idx), 'uint32');
  ret(:) = find(log_idx);
end
