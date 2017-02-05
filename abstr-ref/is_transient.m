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

  elseif strcmp(opt_settings.mode, 'sos')
    % Sos or multi-mode
    res = is_transient_nlin(rec, dyn_list, opt_settings.max_deg);

  elseif strcmp(opt_settings.mode, 'sdsos')
    % Sdsos or multi-mode
    res = is_transient_nlin_sdsos(rec, dyn_list, opt_settings.max_deg);
  end

  res = ;
  % res = ;
