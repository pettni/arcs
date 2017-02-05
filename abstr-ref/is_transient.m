function [res] = is_transient(rec, dyn_list, drec)
  % is_transient: gateway function that delegates transience computation
  % Inputs:
  %  -rec:      hyperrectangle Rec defining a set
  %  -dyn_list: dynamics {f1, f2, f3}
  %  -drec: disturbance

  global opt_settings;

  if length(dyn_list) == 1 && isa(dyn_list{1}, 'cell') && isa(dyn_list{1}{1}, 'double')
    % Linear single-mode
    res = is_transient_lin(rec, dyn_list{1}, drec);
  else
    % Todo: convert linear to poly/yalm
    if strcmp(opt_settings.mode, 'sos')
      % Convert to sdp

      res = is_transient_nlin(rec, dyn_list, drec, opt_settings.max_deg);

    elseif strcmp(opt_settings.mode, 'sdsos')

      new_dyn_list = cell(1, length(dyn_list));
      for m = 1:length(dyn_list)
        fx = dyn_list{m};
        if isa(fx, 'cell') && isa(fx{1}, 'double')
          % Convert to Polynomial
          A_all = [fx{2} flip([fx{1} fx{3}], 2)];  % K A E  in grlex
          fx_new = [];
          for i = 1:size(A_all,1)
            fx_new = [fx_new; Polynomial(size(A_all, 2)-1, A_all(i,:))];
          end
          new_dyn_list{m} = fx_new;
        else
          % Just copy
          new_dyn_list{m} = fx;
        end
      end

      res = is_transient_nlin_sdsos(rec, new_dyn_list, drec, opt_settings.max_deg);
    end 
  end
