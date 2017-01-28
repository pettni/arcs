function result = isTransNLin(rec1, rec2, dyn, deg)
  % isTransLin(rec1, rec2, dyn): Return false if a certificate that 
  % guarantees the non-existence of a flow from rec1 to rec2 under 
  % nonlinear dynamics dyn is found, true otherwise
  % 
  % Inputs 
  %   - rec1, rec2: sets
  %   - dyn = {fx, x, d, drec}: dynamics \dot x = fx(x,d), d \in drec
  %   - deg: relaxation order for moment hierarchy
  global ops;

  isect_rec = intersect(rec1,rec2);
  result = true;

  % overlapping?
  if isect_rec.isFullDim
    warning('isTransLinRec: was called with overlapping Recs - returning true');
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
    h1(flatdim) = 1;
  else
    h1(flatdim) = -1;
  end 

  normal_flow = -h1 * dyn{1};

  if length(dyn) > 2
    % There is disturbance
    F = [isect_rec.xmin' <= dyn{2} <= isect_rec.xmax';
         dyn{4}.xmin' <= dyn{3} <= dyn{4}.xmax'];
  else
    F = [isect_rec.xmin' <= dyn{2} <= isect_rec.xmax'];
  end
  
  [Fnew,obj] = momentmodel(F, normal_flow, deg, 1);
  diagnosis = solvesdp(Fnew, obj, ops);

  if (diagnosis.problem==0 || diagnosis.problem==4) & value(normal_flow)>0
    result = false;  % there is no flow
  end
end