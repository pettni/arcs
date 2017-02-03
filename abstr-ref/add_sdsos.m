function prob = add_sdsos(prob, s1, s2)
  % Add an sdsos constraint to MOSEK 'prob' struct 
  % for the variables x(s1:s2)

  if size(prob.a,1) ~= length(prob.buc)
    error('dimension mismatch')
  end

  [numcon, numvar] = size(prob.a);
  
  L = 1+s2-s1;  % number of variables

  if s1+L > numvar
    error('not enough variables')
  end

  n = (sqrt(1+8*L) - 1)/2;
  if ceil(n) ~= n
    error(strcat('L=', num2str(L), ' does not represent a symmetric matrix'))
  end

  D = sdd_mat(L);

  [d1 d2] = size(D);

  % A = [ -Aold-   0 
  %      O  -I  0  D]
  % with I at pos s1:s1+L

  prob.a = [prob.a            sparse(numcon, d2); 
            sparse(d1, s1-1)  -speye(d1)  sparse(d1, numvar-d1-s1+1)  D];
  prob.buc = [prob.buc; zeros(d1, 1)];
  prob.blc = [prob.blc; zeros(d1, 1)];

  num_cones = d2/3;

  % 1 represents rotated cone 2 x1 x2 >= x3 (var res.symbcon.MSK_CT_RQUAD)
  try
    prob.cones.type = [prob.cones.type ones(1, num_cones)];
    prob.cones.subptr = [prob.cones.subptr length(prob.cones.sub)+(1:3:3*num_cones)];
    prob.cones.sub = [prob.cones.sub numvar+(1:3*num_cones)];
  catch
    prob.cones.type = ones(1, num_cones);
    prob.cones.sub = numvar+(1:3*num_cones);
    prob.cones.subptr = 1:3:3*num_cones;
  end
end
