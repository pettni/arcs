function ret = count_monomials_eq(n,d)
  % Number of monomials in n variables of degree equal to d
  ret = nchoosek(n+d-1, d);
end