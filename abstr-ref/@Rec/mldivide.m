function recvec = mldivide(rec1, rec2)
  % Return the set difference rec1 \ rec2
  % represented as an array of fully dimensional Rec's
  % 
  % Output: Array of Rec representing set difference. 
  %         Non-fully dimensional Rec's are removed.
  %

  if ~intersects(rec1, rec2)
    recvec = rec1;
    return;
  end

  if length(rec1)>1
    recvec = [];
    for rec = rec1
      recvec = [recvec mldivide(rec, rec2)];
    end
    return;
  end

  if length(rec2)>1
    recvec = rec1;
    for rec = rec2
      recvec = mldivide(recvec, rec);
    end
    return;
  end

  recvec = [];
  split_poly = rec1;

  % Sort the dimensions of rec1 from largest to smallest
  % to make parts as uniform as possible in size
  [~, I] = sort(rec1.xmax - rec1.xmin, 2, 'descend');

  for i=I
    X_min = [-Inf*ones(1,rec2.dim); Inf*ones(1,rec2.dim)];
    X_min(2,i) = rec2.xmin(i);
    X_mid = [-Inf*ones(1,rec2.dim); Inf*ones(1,rec2.dim)];
    X_mid(1,i) = rec2.xmin(i);
    X_mid(2,i) = rec2.xmax(i);
    X_max = [-Inf*ones(1,rec2.dim); Inf*ones(1,rec2.dim)];
    X_max(1,i) = rec2.xmax(i);
    recvec = [recvec intersect(Rec(X_min), split_poly) intersect(Rec(X_max), split_poly)];
    split_poly = intersect(Rec(X_mid), split_poly); 
  end

  recvec = removeNonFullDim(recvec);
end

function keep_list = removeNonFullDim(rec)
  keep_list = [];
  for r = rec
    if r.isFullDim
      keep_list = [keep_list r];
    end
  end
end