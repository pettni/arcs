function result = isTransientNLin(rec1,vField, vars, deg)

% Nonlinear systems transience based on barrier certificate type conditions
% N.O., June, 2013
% Caltech
% N.O., Feb, 2016, nonlinear version
% Michigan

[Bx,coefs,mons] = polynomial(vars,deg,1);

epsilon = 10;
vField = vField{1};

F = [];
polys = [];
sos_mult = [];
msos_coefs = [];
% for i =1:length(vars)
%     len = (rec1.xmax(i)-rec1.xmin(i))/2;
%     mid = (rec1.xmax(i)+rec1.xmin(i))/2;
%     polys = [polys; (vars(i)-mid)^2-len^2]; %<=0;
%     [sm sc] = polynomial(vars,deg-2,0);
%     sdisplay(sm)
%     sos_mult = [sos_mult;sm];
%     msos_coefs = [msos_coefs; sc];
%     F = [F; sos(sm)];
% end

for i =1:length(vars)
    %len = (rec1.xmax(i)-rec1.xmin(i))/2;
    %mid = (rec1.xmax(i)+rec1.xmin(i))/2;
    polys = [polys; rec1.xmin(i)-vars(i)]; %<=0;
    polys = [polys; vars(i)-rec1.xmax(i)];
    [sm1 sc1] = polynomial(vars,deg-1,0);
    [sm2 sc2] = polynomial(vars,deg-1,0);
    %sdisplay(sm)
    sos_mult = [sos_mult;sm1;sm2];
    msos_coefs = [msos_coefs; sc1;sc2;];
    F = [F; sos(sm1);sos(sm2)];
end

global ops;

F = [sos(sos_mult'*polys - jacobian(Bx,vars)*vField - epsilon);F];

diagnostics = solvesos(F, [], ops, [coefs;msos_coefs]);
if diagnostics.problem ~= 0 %0 feas, 1 infeas, other something else
    result = 0;
else
    result = 1;
end

end