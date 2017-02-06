function result = is_transient_nlin(rec1, dyn_list, drec, deg)
  % Determine if the modes dyn_list = {fx1, fx2, ...}, where
  % fx = {f, vars} is transient on rec1 using 
  % relaxation order deg

  global ops;

  n_x = length(dyn_list{1}{1});
  vars = dyn_list{1}{2};
  xvars = vars(1:n_x);
  [Bx, coefs] = polynomial(xvars, deg, 1);

  epsilon = 10;

  F = [];
  polys = [];
  sos_mult = [];
  msos_coefs = [];

  if isempty(drec)
    drec = Rec(zeros(2,0));
  end

  all_rec = rec1 * drec;

  for i=1:length(vars)
    polys = [polys; all_rec.xmin(i)-vars(i)]; %<=0;
    polys = [polys; vars(i)-all_rec.xmax(i)];
    [sm1 sc1] = polynomial(vars,deg-1,0);
    [sm2 sc2] = polynomial(vars,deg-1,0);
    sos_mult = [sos_mult; sm1; sm2];
    msos_coefs = [msos_coefs; sc1; sc2];
    F = [F; sos(sm1); sos(sm2)];
  end

  % Add lyapunov constraints
  for i=1:length(dyn_list)
    F = [F; sos(sos_mult'*polys - jacobian(Bx, xvars)*dyn_list{i}{1} - epsilon)];
  end

  diagnostics = solvesos(F, [], ops, [coefs; msos_coefs]);

  if diagnostics.problem ~= 0 %0 feas, 1 infeas, other something else
    result = false;
  else
    result = true;
  end

end