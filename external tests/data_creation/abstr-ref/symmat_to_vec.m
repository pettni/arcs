function [v] = symmat_to_vec(A)
  %symmat_to_vec: obtain vector representation of a symmetric matrix
  if norm(A - A') > 1e-8
    error('A not symmetric')
  end
  n = size(A,1);
  L = n*(n+1)/2;

  v = zeros(1, L);
  for k = 1:L
    [i,j] = symmat_k_ij(k, L);
    v(k) = A(i,j);
  end
end