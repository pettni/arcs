function ret = count_monomials_leq(n,d)
  % Number of monomials in n variables of degree less than or equal to d
  ret = nchoosek(n+d, d);
end
