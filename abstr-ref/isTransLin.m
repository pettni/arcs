function result = isTransLin(rec1, rec2, fx)
  % Function to check if there exists a transition from some point in rec1
  % to some point in rec2 under the _linear_ vector field fx.
  % The computation is done by evaluating the normal flow at corner points.
  % Return true if there exists a flow and false if one can 
  % find a certificate that guarantees the non-existence.

  isect_rec = intersect(rec1,rec2);

  % overlapping
  if isect_rec.isFullDim
    warning('isTransLinRec: was called with overlapping Recs - returning true');
    result = true;
    return;
  end

  % are they adjacent?
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

  result = false;

  n = length(isect_rec.getFullDims); % number of remaining dimensions
  for i = 1:2^n % iterate over vertices
    x_vert = isect_rec.getVertexI(i);
    if length(fx) == 2
      % no disturbance
      if h1*(fx{1}*x_vert' + fx{2}) > 0
        result = true;  % there is a flow
        return;
      end
    else
      % there is disturbance
      for j = 1:2^length(fx{4}.getFullDims)
        d_vert = fx{4}.getVertexI(j);
        if h1*(fx{1}*x_vert' + ...
               fx{2} + ...
               fx{3}*d_vert') > 0
          result = true;  % there is a flow
          return;
        end
      end
    end
  end
end