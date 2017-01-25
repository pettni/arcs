function vert = getVertices(rec)
  % Returns a matrix where each line is the coordinate of a vertex of rec

  if rec.isEmptySet
    vert = [];
    return;
  end

  if rec.dim == 1
    vert = [rec.xmin rec.xmax];
    return;
  end

  if rec.dim == 2
    vert = [rec.xmin(1) rec.xmin(2);
        rec.xmin(1) rec.xmax(2);
        rec.xmax(1) rec.xmax(2);
        rec.xmax(1) rec.xmin(2)];
    return;
  end

  ind_fulldim = rec.getFullDims;
  n_fulldim = length(ind_fulldim);

  vert = repmat(rec.xmin, 2^n_fulldim, 1);
  
  for i = 2:2^n_fulldim
    act_fulldim = find(de2bi(i-1, n_fulldim));
    vert(i,ind_fulldim(act_fulldim)) = rec.xmax(ind_fulldim(act_fulldim));
  end
end