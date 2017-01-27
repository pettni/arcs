function add_mode(part, fx)
	% create_ts: add fxion corresponding to a dynamical mode
	%
	% Linear dynamics:      \dot x = A x + K + E d
	% Polynomial dynamics:  \dot x = f(x, d), f polynomial
	%                                d \in dist_rec
  % Input:
	%   fx = {A, K, E, dist_rec}  in the linear case (last two optional)
	%   fx = {f, vars, dist_rec}, f is an sdpvar vector in variables 'vars'
  %

  act_n = part.ts.add_action();
  part.act_list{act_n} = fx;

  [trans_fun, trans_out_fun, transient_fun] = get_fcns(fx);

  % Figure out transitions
  for i=1:length(part)
    % Neighbor transitions
    adj = part.get_neighbors(i);
    for j=adj
      if trans_fun(part.cell_list(i), part.cell_list(j))
        part.ts.add_transition(i, j, act_n);
      end
    end

    % Out-of-domain
    if trans_out_fun(part.cell_list(i), part.domain)
      part.ts.add_transition(i, length(part)+1, act_n);
    end

    % Self transitions
    if ~transient_fun(part.cell_list(i))
      part.ts.add_transition(i, i, act_n);
    end
  end

  % Figure out progress groups
  if transient_fun(part.domain)
    % Whole domain progress groups
    part.ts.add_progress_group([act_n], 1:length(part))
  else
    % Individual state progress group
    for i = 1:length(part)
      if transient_fun(part.cell_list(i))
        part.ts.add_progress_group([act_n], [i])
      end
    end
  end
end