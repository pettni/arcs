%% get_fcns: get appropriate function for given dynamics
function [trans_fun, trans_out_fun, transient_fun] = get_fcns(fx)

  if strcmp(class(fx{1}), 'sdpvar')
    dyn_type = 'polynomial';
  else
    dyn_type = 'linear';
  end

  if strcmp(dyn_type, 'linear')
    trans_fun = @(p1, p2) isTransLin(p1, p2, fx);
    trans_out_fun = @(p1, dom) isTransOutLin(p1, dom, fx);
    transient_fun = @(p1) isTransientLin(p1, fx);
  elseif strcmp(dyn_type, 'polynomial')
    % Degrees seem arbitrary
    trans_fun = @(p1, p2) isTransNLin(p1, p2, fx, 2);
    trans_out_fun = @(p1, dom) isTransOutNLin(p1, dom, fx, 2);
    transient_fun = @(p1) isTransientNLin(p1, fx, 4);
  else
    error('unknown dynamics')
  end
