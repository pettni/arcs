function [ret] = sdd_mat(L)
  n = (sqrt(1+8*L) - 1)/2;

  if (n ~= (sqrt(1+8*L) - 1)/2)
    error('L not the length of a symmetric matrix representation')
  end

  num_vars = n*(n-1)/2;

  ret = sparse(L, 3*num_vars);

  for k = 1:L
    [i,j] = symmat_k_ij(k, L);
    if i==j
      for l = 1:n
        if l == i
          continue
        end
        mat_idx = symmat_ij_k(min(i, l), max(i, l)-1, num_vars);
        el_idx =  1 + (i>l);
        ret(k, 3*(mat_idx-1) + el_idx) = 1 + (el_idx == 1);   
                                           % msk cone: 2 x1 x2 >= x3^2
      end
    else
      mat_idx = symmat_ij_k(min(i,j), max(i,j)-1, num_vars);
      el_idx = 3;
      ret(k, 3*(mat_idx-1) + el_idx) = 1;
    end
  end
end
