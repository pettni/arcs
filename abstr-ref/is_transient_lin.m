function result = is_transient_lin(rec1, fx, drec)
  % Determine if the mode fx = {A, K, E} is transient on rec1
  
  global ops;

  Ad = fx{1};
  Kd = fx{2};

  n = size(Ad, 2);

  if length(fx) == 3 && nargin == 3 && ~isempty(drec)
    Ed = fx{3};
  else
    % No disturbance
    Ed = zeros(n,0);
    drec = Rec(zeros(2,0));
  end

  n_d = size(Ed, 2);
  assert (drec.dim == n_d);

  r1 = rank(Ad);

  if r1 == n
    % Full-rank system
    if n_d
      % disturbance
      result = true;
      for j = 1:2^length(drec.getFullDims)
        % for each disturbance vertex
        d_vert = drec.getVertexI(j);
        xe = -Ad\(Kd+Ed*d_vert');
        if isInside(rec1, xe)
          result = false;
          return;
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
    if ~n_d
      % no disturbance
      if rank([Ad Kd]) ~= r1
        result = true;
        return
      end
      Consts = [Ad*x == -Kd, rec1.xmin' <= x, x <= rec1.xmax'];
    else
      d = sdpvar(n_d, 1);
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