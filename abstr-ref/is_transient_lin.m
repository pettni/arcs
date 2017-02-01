function result = is_transient_lin(rec1, fx)
  % Determine if the mode fx = {A, K, (E, drec)} is transient on rec1
  
  global ops;

  Ad = fx{1};
  Kd = fx{2};

  n = size(Ad, 2);

  if length(fx) == 4
    % there is disturbance
    Ed = fx{3};
    n_d = size(Ed, 2);
  end

  r1 = rank(Ad);

  if r1 == n
    % Full-rank system
    if length(fx) == 4
      % disturbance
      result = true;
      for j = 1:2^length(fx{4}.getFullDims)
        % for each disturbance vertex
        d_vert = fx{4}.getVertexI(j);
        xe = -Ad\(Kd+Ed*d_vert');
        if isInside(rec1, xe)
          result = false;
          break;
        end
      end
    else
      % no disturbance
      xe = -Ad\Kd;
      result = not(isInside(rec1, xe));
    end
  else
    x = sdpvar(n,1);
    Obj = ones(1,n)*x;
    if length(fx) == 2
      % no disturbance
      if rank([Ad Kd]) ~= r1
        result = true;
        return
      end
      Consts = [Ad*x == -Kd, rec1.xmin' <= x, x <= rec1.xmax'];
    else
      Ed = fx{3};
      n_d = size(Ed, 2);
      d = sdpvar(n_d, 1)
      drec = fx{4};
      Consts = [Ad*x + Ed * d == -Kd, rec1.xmin' <= x, x <= rec1.xmax', ...
                                      drec.xmin' <= d <= drec.xmax'];
    end
    diagnostics = solvesdp(Consts, Obj, ops);
    if diagnostics.problem ~= 0 %0 feas, 1 infeas, other something else
      result = true;
    else
      result = false;
    end
  end