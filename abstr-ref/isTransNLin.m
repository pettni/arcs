function result = isTransNLin(rec1, rec2, fx, deg)

  % Function to check if there exists a transition from some point in rec1
  % to some point in rec2 under the _linear_ vector field fx.
  % The computation is done by evaluating the normal flow at corner points.
  % Return true if there exists a flow and false if one can 
  % find a certificate that guarantees the non-existence.
  global ops;

  isect_rec = intersect(rec1,rec2);
  result = true;

  % overlapping?
  if isect_rec.isFullDim
    warning('isTransLinRec: was called with overlapping Recs - returning true');
    result = true;
    return;
  end

  % adjacent?
  if length(isect_rec.getFlatDims) ~= 1
    result = false;
    return;
  end

  % find "flat" dimension
  flatdim = isect_rec.getFlatDims;
  h1 = zeros(1,rec1.dim);
  if rec1.xmax(flatdim)==rec2.xmin(flatdim)
    h1(flatdim)=1;
  else
    h1(flatdim)=-1;
  end 
  % h1 normal vector

  poly = -h1*(fx{1});

  F = isect_rec.xmin' <= fx{2} <= isect_rec.xmax';
  
  [Fnew,obj] = momentmodel(F,poly,deg,1);
  diagnosis = solvesdp(Fnew, obj, ops);

  if (diagnosis.problem==0 || diagnosis.problem==4) & value(poly)>0
    result = false;  % there is no flow
  end
end