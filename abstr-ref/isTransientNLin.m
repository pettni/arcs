function result = isTransientNLin(rec1, act, deg)

% Nonlinear systems transience based on barrier certificate type conditions
% N.O., June, 2013
% Caltech
% N.O., Feb, 2016, nonlinear version
% Michigan

global ops;

fx = act{1};
vars = act{2};

[Bx,coefs,mons] = polynomial(vars,deg,1);

epsilon = 10;
fx = fx{1};

F = [];
polys = [];
sos_mult = [];
msos_coefs = [];

for i =1:length(vars)
    polys = [polys; rec1.xmin(i)-vars(i)]; %<=0;
    polys = [polys; vars(i)-rec1.xmax(i)];
    [sm1 sc1] = polynomial(vars,deg-1,0);
    [sm2 sc2] = polynomial(vars,deg-1,0);
    sos_mult = [sos_mult;sm1;sm2];
    msos_coefs = [msos_coefs; sc1;sc2;];
    F = [F; sos(sm1); sos(sm2)];
end

F = [F; sos(sos_mult'*polys - jacobian(Bx,vars)*fx - epsilon)];

diagnostics = solvesos(F, [], ops, [coefs; msos_coefs]);
if diagnostics.problem ~= 0 %0 feas, 1 infeas, other something else
    result = 0;
else
    result = 1;
end

end