function abstract(part, dyn_list, vars, drec)
    % abstract: abstract the dynamics in dyn_list = {fx1, ..., fxM}
    %
    % Dynamics:  \dot x = f_m(x, d), f polynomial
    %                                d \in drec
    % Input:
    %   fx : sdpvar vector in variables 'vars'
    %   vars : sdpvars
    %   drec : bound on disturbance variables (must be at the end of vars)

  global opt_settings
  
  if nargin < 4
    drec = [];
  end

  M = length(dyn_list);

  part.create_ts();
  part.d_rec = drec;

  % Convert and save in partition
  for m = 1:M
    act_num = part.ts.add_action();
    fx = dyn_list{m};
    % Linear?
    if degree(fx, vars) <= 1
      n_x = length(fx);
      A_all = zeros(n_x, 1 + length(vars));
      for i = 1:n_x
        [coef, mons] = coefficients(fx(i), vars);
        for j = 1:length(coef)
          b = mono_rank_grlex(length(vars), degree(mons(j), flip(vars)));
          A_all(i, b) = coef(j);
        end
      end
      K = A_all(:, 1);
      A = A_all(:, 2:1+n_x);
      E = A_all(:, 2+n_x:end);
      part.dyn_list{act_num} = {A, K, E};
    elseif strcmp(opt_settings.mode, 'sdsos')
      % Convert to Polynomial
      part.dyn_list{act_num} = Polynomial(fx, vars);
    elseif strcmp(opt_settings.mode, 'sos')
      % Add as-is
      part.dyn_list{act_num} = {fx, vars};
    else
      error('opt_settings.mode must be sdsos or sos')
    end
    part.dyn_list_orig{act_num} = {dyn_list{m}, vars};
  end

  % Figure out transitions
  for m = 1:M
    for i=1:length(part)
      % Neighbor transitions
      adj = part.get_neighbors(i);
      for j=adj
        if is_trans(part.cell_list(i), part.cell_list(j), part.dyn_list{m}, drec)
          part.ts.add_transition(i, j, m);
        end
      end

      % Out-of-domain
      if is_trans_out(part.cell_list(i), part.domain, part.dyn_list{m}, drec)
        part.ts.add_transition(i, length(part)+1, m);
      end

      % Self transitions
      if ~is_transient(part.cell_list(i), {part.dyn_list{m}}, drec)
        part.ts.add_transition(i, i, m);
      end
    end
  end

end