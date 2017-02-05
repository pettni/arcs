function [res] = is_trans(rec1, rec2, dyn)
  % is_trans: determine if there is a transition from rec1 to rec2
  % under mode dyn

  
  if ~strcmp(class(dyn{2}), 'sdpvar')
    % Linear
    res = is_trans_lin(rec1, rec2, dyn);
  else
    % res = is_trans_nlin(rec1, rec2, dyn, 4);
    res = is_trans_nlin_sdsos(rec1, rec2, dyn, 4);
  end
