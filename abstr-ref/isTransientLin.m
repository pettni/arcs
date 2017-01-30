function result = isTransientLin(rec1, fx)

% For linear systems transience is equivalent to nonexistence of eq. pnts.
% under the assumption that the system does no involve modes on the unit
% circle
% N.O., June, 2013
% Caltech

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
    r2 = rank([Ad Kd]);
    if r2~=r1
      result = true;
      return
    else
      Hr = [eye(rec1.dim); -eye(rec1.dim)]; Kr = [rec1.xmax'; rec1.xmin'];
      x = sdpvar(n,1);
      Consts = [Ad*x == -Kd, Hr*x <= Kr];
      Obj = ones(1,n)*x;
      global ops;
      diagnostics = solvesdp(Consts, Obj, ops);
      if diagnostics.problem ~= 0 %0 feas, 1 infeas, other something else
        result = true;
      else
        result = false;
      end
    end
end