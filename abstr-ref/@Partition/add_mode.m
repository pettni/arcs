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

  % Figure out transitions
  for i=1:length(part)
    % Neighbor transitions
    adj = part.get_neighbors(i);
    for j=adj
      if is_trans(part.cell_list(i), part.cell_list(j), fx)
        part.ts.add_transition(i, j, act_n);
      end
    end

    % Out-of-domain
    if is_trans_out(part.cell_list(i), part.domain, fx)
      part.ts.add_transition(i, length(part)+1, act_n);
    end

    % Self transitions
    if ~is_transient(part.cell_list(i), {fx})
      part.ts.add_transition(i, i, act_n);
    end
  end

end