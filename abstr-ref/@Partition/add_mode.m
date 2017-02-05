function add_mode(part, fx, vars, drec)
	% create_ts: add action corresponding to a dynamical mode 'fx'
	%
	% Polynomial dynamics:  \dot x = f(x, d), f polynomial
	%                                d \in drec
  % Input:
	%   fx : sdpvar vector in variables 'vars'
  %   drec : bound on disturbance variables (must be at the end of vars)
  %

  global opt_settings

  if nargin < 4
    drec = [];
  end

  act_n = part.ts.add_action();

  if strcmp(opt_settings.mode, 'sdsos')
    % Convert to Polynomial
    part.dyn_list{act_n} = {Polynomial(fx, vars), drec};
  elseif strcmp(opt_settings.mode, 'sos')
    part.dyn_list{act_n} = {fx, vars, drec};
  else
    error('opt_settings.mode must be sdsos or sos')
  end

  % Figure out transitions
  for i=1:length(part)
    % Neighbor transitions
    adj = part.get_neighbors(i);
    for j=adj
      if is_trans(part.cell_list(i), part.cell_list(j), part.dyn_list{act_n})
        part.ts.add_transition(i, j, act_n);
      end
    end

    % Out-of-domain
    if is_trans_out(part.cell_list(i), part.domain, part.dyn_list{act_n})
      part.ts.add_transition(i, length(part)+1, act_n);
    end

    % Self transitions
    if ~is_transient(part.cell_list(i), {part.dyn_list{act_n}})
      part.ts.add_transition(i, i, act_n);
    end
  end

end