function result = isTransientNLinD(rec1,vField, vars, deg, dbound)

% Nonlinear systems transience based on barrier certificate type conditions
% N.O., June, 2013
% Caltech
% N.O., Feb, 2016, nonlinear version
% Michigan

vField = vField{1};
dim = size(vField,1);

varsx = vars(1:dim);
varsd = vars(dim+1:end);

[Bx,coefs,mons] = polynomial(varsx,deg,1);

epsilon = 10;

F = [];
polys = [];
sos_mult = [];
msos_coefs = [];


for i =1:length(varsx)
    %len = (rec1.xmax(i)-rec1.xmin(i))/2;
    %mid = (rec1.xmax(i)+rec1.xmin(i))/2;
    polys = [polys; rec1.xmin(i)-varsx(i)]; %<=0;
    polys = [polys; varsx(i)-rec1.xmax(i)];
    [sm1 sc1] = polynomial(vars,deg-1,0);
    [sm2 sc2] = polynomial(vars,deg-1,0);
    %sdisplay(sm)
    sos_mult = [sos_mult;sm1;sm2];
    msos_coefs = [msos_coefs; sc1;sc2;];
    F = [F; sos(sm1);sos(sm2)];
end

for i =1:length(varsd)
    polys = [polys; -dbound-varsd(i)]; %<=0;
    polys = [polys; varsx(i)-dbound];
    [sm1 sc1] = polynomial(vars,deg-1,0);
    [sm2 sc2] = polynomial(vars,deg-1,0);
    %sdisplay(sm)
    sos_mult = [sos_mult;sm1;sm2];
    msos_coefs = [msos_coefs; sc1;sc2;];
    F = [F; sos(sm1);sos(sm2)];
end

global ops;

F = [sos(sos_mult'*polys - jacobian(Bx,varsx)*vField - epsilon);F];

diagnostics = solvesos(F, [], ops, [coefs;msos_coefs]);
if diagnostics.problem ~= 0 %0 feas, 1 infeas, other something else
    result = 0;
else
    result = 1;
end

end