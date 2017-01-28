function result = isTransientNLin(rec1, dyn_list, deg)
  % Nonlinear systems transience based on barrier certificate type conditions
  % N.O., June, 2013
  % Caltech
  % N.O., Feb, 2016, nonlinear version
  % Michigan

  global ops;

  [Bx, coefs, mons] = polynomial(dyn_list{1}{2}, deg, 1);

  epsilon = 10;

  F = [];
  polys = [];
  sos_mult = [];
  msos_coefs = [];

  for dyn_nr=1:length(dyn_list)
    dyn = dyn_list{dyn_nr};
    if length(dyn) > 2
      % Disturbance
      all_rec = Rec([rec1.xmin dyn{4}.xmin; rec1.xmax dyn{4}.xmax]);
      all_vars = [dyn{2}; dyn{3}];
    else
      % No disturbance
      all_rec = rec1;
      all_vars = dyn{2};
    end

    % Add bounding conditions
    for i =1:length(all_vars)
        polys = [polys; all_rec.xmin(i)-all_vars(i)]; %<=0;
        polys = [polys; all_vars(i)-all_rec.xmax(i)];
        [sm1 sc1] = polynomial(all_vars,deg-1,0);
        [sm2 sc2] = polynomial(all_vars,deg-1,0);
        sos_mult = [sos_mult;sm1;sm2];
        msos_coefs = [msos_coefs; sc1;sc2;];
        F = [F; sos(sm1); sos(sm2)];
    end

    F = [F; sos(sos_mult'*polys - jacobian(Bx, dyn{2})*dyn{1} - epsilon)];
  end

  diagnostics = solvesos(F, [], ops, [coefs; msos_coefs]);

  if diagnostics.problem ~= 0 %0 feas, 1 infeas, other something else
      result = false;
  else
      result = true;
  end

end