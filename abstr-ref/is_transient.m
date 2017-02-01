function [res] = is_transient(rec, dyn_list)
  % is_transient: gateway function that delegates transience computation
  % Inputs:
  %  -rec:      hyperrectangle Rec defining a set
  %  -dyn_list: dynamics {f1, f2, f3}, where fi = {A, B, (E, d_rec)} or
  %                                          fi = {f, xvar, dvar, drec}
  % 

  if length(dyn_list) == 1 && ~strcmp(class(dyn_list{1}{2}), 'sdpvar')
    % Linear single-mode
    res = is_transient_lin(rec, dyn_list{1});
    return
  end

  % Find an x variable
  xvar = [];
  for i = 1:length(dyn_list)
    if strcmp(class(dyn_list{i}{2}), 'sdpvar')
      xvar = dyn_list{i}{2};
      break;
    end
  end
  if isempty(xvar)
    % All are linear, create new variable
    xvar = sdpvar(size(dyn_list{1}{1}, 2), 1);
  end

  % We have multi-mode system: convert linear ones
  new_dyn_list = {};
  for i = 1:length(dyn_list)
    fx = dyn_list{i};
    if ~strcmp(class(fx{2}), 'sdpvar')
      % Linear
      A = fx{1};
      K = fx{2};
      if length(fx) == 2
        % No disturbance
        new_dyn_list{i} = {A * xvar + K, xvar};
      else
        E = fx{3};
        d_rec = fx{4};
        dvar = sdpvar(size(E,2), 1);
        new_dyn_list{i} = {A * xvar + K + E * dvar, xvar, dvar, d_rec};
      end
    else
      % Nonlinear--just copy
      new_dyn_list{i} = fx;
    end
  end
  res = is_transient_nlin(rec, new_dyn_list, 4);
