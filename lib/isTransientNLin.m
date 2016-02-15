function result = isTransientNLin(rec1,vField, vars, deg)

% Nonlinear systems transience based on barrier certificate type conditions
% N.O., June, 2013
% Caltech
% N.O., Feb, 2016, nonlinear version
% Michigan

[Bx,coefs,mons] = polynomial(vars,deg,0);

epsilon = 0.1;
vField = vField{1};

F = [];
polys = [];
sos_mult = [];
sos_coefs = [];
for i =1:length(vars)
    len = (rec1.xmax(i)-rec1.xmin(i))/2;
    mid = (rec1.xmax(i)+rec1.xmin(i))/2;
    polys = [polys; (vars(i)-mid)^2-len^2]; %<=0;
    [sm sc] = polynomial(vars,deg-2,0);
    sos_mult = [sos_mult;sm];
    sos_coefs = [sos_coefs; sc];
    F = [F; sos(sos_mult)];
end

F = [sos(sos_mult'*polys - jacobian(Bx,vars)*vField - eps);F];

diagnostics = solvesos(F, [coefs;sos_coefs]);
if diagnostics.problem ~= 0 %0 feas, 1 infeas, other something else
    result = 0;
else
    result = 1;
end

end