function [res] = is_trans(rec1, rec2, dyn)
  % is_trans: determine if there is a transition from rec1 to rec2
  % under mode dyn

  global opt_mode
  
  if ~strcmp(class(dyn{2}), 'sdpvar')
    % Linear
    res = is_trans_lin(rec1, rec2, dyn);
  else
    if strcmp(opt_settings.mode, 'sos')
      res = is_trans_nlin(rec1, rec2, dyn, opt_settings.max_deg);
    elseif strcmp(opt_settings.mode, 'sdsos')
      res = is_trans_nlin_sdsos(rec1, rec2, dyn, opt_settings.max_deg);
    end
  end
