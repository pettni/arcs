function [res] = is_trans(rec1, rec2, dyn, drec)
  % is_trans: determine if there is a transition from rec1 to rec2
  % under mode dyn = {A, K, E}  (linear)
  %                = {fx, vars} (yalmip-sos)
  %                = pol        (Polynomial-sdsos)
  global opt_settings
  
  if isa(dyn, 'cell') && isa(dyn{1}, 'double')
    res = is_trans_lin(rec1, rec2, dyn, drec);
  elseif isa(dyn, 'cell') && isa(dyn{1}, 'sdpvar')
    res = is_trans_nlin(rec1, rec2, dyn, drec, opt_settings.max_deg);
  elseif isa(dyn, 'Polynomial')
    res = is_trans_nlin_sdsos(rec1, rec2, dyn, drec, opt_settings.max_deg);
  end
