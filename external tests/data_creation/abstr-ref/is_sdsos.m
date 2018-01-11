function [ret, M_list] = is_sdsos(A)
  % is_sdsos: determine if a matrix A is sdsos
  % A can be a square matrix or the (n*(n+1)/2) row vector representation of
  % a symmetric matrix

  % Set up a conelp and see if it is feasible
  % [I] [c] = [v];  c sdsos

  if size(A,1) == 1
    v = A;
  else
    if norm(A-A') > 1e-10
      error('A is not symmetric')
    end
    n = size(A, 1);
    v = zeros(1, n*(n+1)/2);

    for i = 1:length(v)
      [j, k] = symmat_k_ij(i, length(v));
      v(i) = A(j,k);
    end
  end
  n = (sqrt(1+8*length(v)) - 1)/2;
  if floor(n) ~= n
    error('vector does not represent a symmetric matrix')
  end
  
  prob = [];
  prob.a = speye(length(v));
  prob.buc = v';
  prob.blc = v';

  prob = add_sdsos(prob, 1, length(v));

  prob.c = zeros(size(prob.a, 2), 1);
  prob.c(length(v)+1:end) = 1;

  [~, res] = mosekopt('minimize echo(0)', prob);
  M = res.sol.itr.xx(length(v)+1:end);

  M_list = cell(1, length(M)/3);
  for i=1:3:length(M)
    Mlift = zeros(n);
    k = 1+(i-1)/3;
    [id, jd] = symmat_k_ij(k, n*(n-1)/2);
    Mlift(id, id) = 2*M(i);  % Mosek cone compensation
    Mlift(id, jd+1) = M(i+2);
    Mlift(jd+1, id) = M(i+2);
    Mlift(jd+1, jd+1) = M(i+1);
    if det(Mlift) < 0 || trace(Mlift) < 0
      error('not pos def')
    end
    M_list{k} = Mlift;
  end

  if strcmp(res.sol.itr.prosta, 'PRIMAL_AND_DUAL_FEASIBLE')
    ret = true;
  else
    ret = false;
  end
end