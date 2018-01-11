function [i, j] = symmat_k_ij(k, L)
  % symmat_ij_k(i, j, L): Given a symmetric matrix Q represented by a vector
  % V = [Q_11, ... Q_1n, Q_22, ..., ... Q_nn] of length L,
  % for given k compute i,j s.t. Q_ij = V(k)
  n = (sqrt(1+8*L) - 1)/2;

  if (n ~= (sqrt(1+8*L) - 1)/2)
    error('L not the length of a symmetric matrix representation')
  end

  % # get first index
  i = ceil( (2*n+1)/2 - sqrt( ((2*n+1)/2)^2 - 2*k ));

  % # second index
  k1 = (2*n+2-i)*(i-1)/2 - 1;
  j = i-1 + k - k1 - 1;
end
