function result = isTransientLinMulti(rec1, dyn_list, deg)
  % isTransientLinMulti: Multi-field transience for linear system
  % redirects to appropriate function
  if length(dyn_list) == 1
    % Single mode: exact
    result = isTransientLin(rec1, dyn_list{1});
  else
    % Multi-mode
    new_dyn_list = {};
    if length(dyn_list{1}) == 2
      % No disturbance
      n_x = size(dyn_list{1}{1}, 2);
      x = sdpvar(n_x, 1);
      for i=1:length(dyn_list)
        new_dyn_list{i} = {dyn_list{i}{1}*x + dyn_list{i}{2}, x};
      end
    else
      error('not implemented')
    end
    result = isTransientNLinMulti(rec1, new_dyn_list, deg);
  end

end