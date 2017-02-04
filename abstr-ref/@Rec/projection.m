function prec = projection(rec, dims)
  % Returns the projection of rec onto dimensions dims
  if max(dims) > rec.dim
    error('Dimension too large')
  end
  pxmin = rec.xmin(dims);
  pxmax = rec.xmax(dims);
  prec = Rec([pxmin; pxmax]);
end