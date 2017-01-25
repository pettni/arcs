function frec = getFacetI(rec,i)
  % return facet number i of rec
  if i<1 || i>2*rec.dim
    error('index out of range')
  end
  fmin = rec.xmin;
  fmax = rec.xmax;
  chind = 1+mod(i-1, rec.dim)
  if i <= rec.dim
    fmax(chind) = fmin(chind);
  else
    fmin(chind) = fmax(chind);
  end
  frec = Rec([fmin; fmax]);
end