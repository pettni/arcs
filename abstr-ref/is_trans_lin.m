function result = is_trans_lin(rec1, rec2, dyn, drec)
  % is_trans_lin(rec1, rec2, dyn, drec): Return false if a certificate that 
  % guarantees the non-existence of a flow from rec1 to rec2 under 
  % linear dynamics dyn is found, true otherwise
  % 
  % Inputs:
  %   - rec1, rec2: sets
  %   - dyn = {A, K, E}: dynamics \dot x = Ax + K + Ed, d \in drec
  %   - drec: disturbance Rec

  isect_rec = intersect(rec1,rec2);

  if nargin < 4 || isempty(drec)
    drec = Rec(zeros(2,0));
  end

  % overlapping
  if isect_rec.isFullDim
    warning('is_trans_lin: was called with overlapping Recs - returning true');
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
    if length(dyn) == 2
      % no disturbance
      if h1*(dyn{1}*x_vert' + dyn{2}) > 0
        result = true;  % there is a flow
        return;
      end
    else
      % there is disturbance
      for j = 1:2^length(drec.getFullDims)
        d_vert = drec.getVertexI(j);
        if h1*(dyn{1}*x_vert' + ...
               dyn{2} + ...
               dyn{3}*d_vert') > 0
          result = true;  % there is a flow
          return;
        end
      end
    end
  end
end