function [ret] = symmat_ij_k(i, j, L)
  % symmat_ij_k(i, j, L): Given a symmetric matrix Q represented by a vector
  % V = [Q_11, ... Q_1n, Q_22, ..., ... Q_nn] of length L,
  % for given i,j , compute k s.t. Q_ij = V(k)

  n = (sqrt(1+8*L) - 1)/2;  % side of matrix

  if (n ~= (sqrt(1+8*L) - 1)/2)
    error('L not the length of a symmetric matrix representation')
  end

  i_at1 = min(i, j);
  j_at1 = max(i, j);
  ret = (n + n-i_at1)*(i_at1-1)/2 + j_at1;
