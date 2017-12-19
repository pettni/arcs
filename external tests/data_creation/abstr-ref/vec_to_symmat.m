%% vec_to_symmat: convert a vector representation of a symmetric matrix
function [AA] = vec_to_symmat(v)
  n = (sqrt(1+8*length(v)) - 1)/2;
  AA = zeros(n, n);
  for k = 1:length(v)
    [i,j] = symmat_k_ij(k, length(v));
    AA(i,j) = v(k);
    AA(j,i) = v(k);
  end
end