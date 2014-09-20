function result = isTransientLin(rec1,vField)

% For linear systems transience is equivalent to nonexistence of eq. pnts.
% under the assumption that the system does no involve modes on the unit
% circle
% N.O., June, 2013
% Caltech

Ad = vField.A;
Kd = vField.K;

n = length(Ad);

r1 = rank(Ad);

if r1 == n
    xe = -Ad\Kd;
    result = not(isInside(rec1,xe));
    return
else
    r2 = rank([Ad Kd]);
    if r2~=r1
        result = 1;
        return
    else
        Hr = [eye(rec1.dim); -eye(rec1.dim)]; Kr = [rec1.xmax'; rec1.xmin'];
        x = sdpvar(n,1);
        Consts = [Ad*x == -Kd, Hr*x <= Kr];
        Obj = ones(1,n)*x;
        global ops;
        diagnostics = solvesdp(Consts, Obj, ops);
        if diagnostics.problem ~= 0 %0 feas, 1 infeas, other something else
            result = 1;
        else
            result = 0;
        end
    end
end