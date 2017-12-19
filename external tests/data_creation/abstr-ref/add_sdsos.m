function prob = add_sdsos(prob, s1, s2, eps)
  % Add an sdsos constraint to MOSEK 'prob' struct 
  % for the variables x(s1:s2)
  % if eps specified, add additional margin A - I is sdsos

  if size(prob.a,1) ~= length(prob.buc)
    error('dimension mismatch')
  end

  if nargin < 4
    eps = 0;
  end

  [numcon, numvar] = size(prob.a);
  
  L = 1+s2-s1;  % number of variables
  n = (sqrt(1+8*L) - 1)/2;

  if s2 > numvar
    error('not enough variables')
  end

  n = (sqrt(1+8*L) - 1)/2;
  if ceil(n) ~= n
    error(strcat('L=', num2str(L), ' does not represent a symmetric matrix'))
  end

  if ~isfield(prob, 'cones')
    prob.cones = [];
  end

  if ~isfield(prob.cones, 'type')
    prob.cones.type = [];
    prob.cones.sub = [];
    prob.cones.subptr = [];
  end

  num_vars = n*(n-1)/2;
  D = sparse(L, 3*num_vars);
  for k = 1:L
    [i,j] = symmat_k_ij(k, L);
    if i==j
      for l = 1:n
        if l == i
          continue
        end
        mat_idx = symmat_ij_k(min(i, l), max(i, l)-1, num_vars);
        el_idx =  1 + (i>l);
        D(k, 3*(mat_idx-1) + el_idx) = 1 + (el_idx == 1);   
                                           % msk cone: 2 x1 x2 >= x3^2
      end
    else
      mat_idx = symmat_ij_k(min(i,j), max(i,j)-1, num_vars);
      el_idx = 3;
      D(k, 3*(mat_idx-1) + el_idx) = 1;
    end
  end

  [d1 d2] = size(D);

  % A = [ -Aold-   0 
  %      O  -I  0  D]
  % with I at pos s1:s1+L

  prob.a = [prob.a            sparse(numcon, d2); 
            sparse(d1, s1-1)  -speye(d1)  sparse(d1, numvar-d1-s1+1)  D];
  prob.buc = [prob.buc; -eps*symmat_to_vec(eye(n))'];
  prob.blc = [prob.blc; -eps*symmat_to_vec(eye(n))'];

  num_cones = d2/3;

  % 1 represents rotated cone 2 x1 x2 >= x3 (var res.symbcon.MSK_CT_RQUAD)
  prob.cones.type = [prob.cones.type ones(1, num_cones)];
  prob.cones.subptr = [prob.cones.subptr length(prob.cones.sub)+(1:3:3*num_cones)];
  prob.cones.sub = [prob.cones.sub numvar+(1:3*num_cones)];
end
