function vert = getVertexI(rec, i)
  % Return vertex number i of rec; good for iterating over vertices without using memory
  if i == 1
    vert = rec.xmin;
    return;
  end

  ind_fulldim = rec.getFullDims;
  n_fulldim = length(ind_fulldim);
  if (i > 2^n_fulldim)
    error('index too high')
  end

  vert = rec.xmin;
  
  act_fulldim = logical(dec2bin(i-1, n_fulldim) - '0');
  vert(ind_fulldim(act_fulldim)) = rec.xmax(ind_fulldim(act_fulldim));
  
end