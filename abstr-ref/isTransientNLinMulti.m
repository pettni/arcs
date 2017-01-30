function result = isTransientNLin(rec1, dyn_list, deg)
  % Nonlinear systems transience based on barrier certificate type conditions
  % N.O., June, 2013
  % Caltech
  % N.O., Feb, 2016, nonlinear version
  % Michigan

  global ops;

  [Bx, coefs] = polynomial(dyn_list{1}{2}, deg, 1);

  epsilon = 10;

  F = [];
  polys = [];
  sos_mult = [];
  msos_coefs = [];

  xvars = dyn_list{1}{2};
  all_vars = xvars;
  for i=1:length(dyn_list)
    if length(dyn_list{i}) > 2
      all_vars = [all_vars; dyn_list{i}{3}];
    end
  end
  
  % Add x constraints
  for i=1:length(xvars)
    polys = [polys; rec1.xmin(i)-xvars(i)]; %<=0;
    polys = [polys; xvars(i)-rec1.xmax(i)];
    [sm1 sc1] = polynomial(all_vars,deg-1,0);
    [sm2 sc2] = polynomial(all_vars,deg-1,0);
    sos_mult = [sos_mult; sm1; sm2];
    msos_coefs = [msos_coefs; sc1; sc2];
    F = [F; sos(sm1); sos(sm2)];
  end

  % Add d constraints
  for i=1:length(dyn_list)
    if length(dyn_list{i}) > 2
      d_vars = dyn_list{3};
      d_rec = dyn_list{4};
      for j=1:length(d_var)
        polys = [polys; d_rec.xmin(i)-d_vars(i)]; %<=0;
        polys = [polys; d_vars(i)-d_rec.xmax(i)];
      end
      [sm1 sc1] = polynomial(all_vars,deg-1,0);
      [sm2 sc2] = polynomial(all_vars,deg-1,0);
      sos_mult = [sos_mult; sm1; sm2];
      msos_coefs = [msos_coefs; sc1; sc2];
      F = [F; sos(sm1); sos(sm2)];
    end
  end

  % Add lyapunov constraints
  for i=1:length(dyn_list)
    F = [F; sos(sos_mult'*polys - jacobian(Bx, xvars)*dyn_list{i}{1} - epsilon)];
  end

  diagnostics = solvesos(F, [], ops, [coefs; msos_coefs]);
  % value(coefs)

  if diagnostics.problem ~= 0 %0 feas, 1 infeas, other something else
    result = false;
  else
    result = true;
  end

end