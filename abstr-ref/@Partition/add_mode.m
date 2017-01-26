function add_mode(part, act, dist_level)
	% create_ts: add action corresponding to a dynamical mode
	%
	% Linear dynamics:      \dot x = A x + K + E d
	% Polynomial dynamics:  \dot x = f(x, d), f polynomial
	% 
  % Input:
	%   act = {A, K, E}  in the linear case
	%   act = {f, vars}, f is an sdpvar vector in variables 'vars'
  %

	if strcmp(class(act{1}), 'sdpvar')
    dyn_type = 'polynomial';
  else
    dyn_type = 'linear';
  end

  disturbance = nargin == 3;

  if strcmp(dyn_type, 'linear') && disturbance
    error('Linear case with disturbance not implemented')
  end

  act_n = part.ts.add_action();
  part.act_list{act_n} = {act, dyn_type, disturbance};  % save dyn info

  if strcmp(dyn_type, 'linear') && ~disturbance
    trans_fun = @(p1, p2) isTransLin(p1, p2, act);
    trans_out_fun = @(p1) isTransOutLin(p1, part.domain, act);
    transient_fun = @(p1) isTransientLin(p1, act);
  end

  if strcmp(dyn_type, 'polynomial') && ~disturbance
    % Degrees seem arbitrary
    trans_fun = @(p1, p2) isTransNLin(p1, p2, act, 2);
    trans_out_fun = @(p1) isTransOutNLin(p1, part.domain, act, 2);
    transient_fun = @(p1) isTransientNLin(p1, act, 4);
  end

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
    if trans_out_fun(part.cell_list(i))
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