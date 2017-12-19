function [ret, dim] = isNeighbor(rec1, rec2)
  % Returns true if rec1 and rec2 are neighbors such that their intersection
  % loses exactly one dimension.
  % Outputs:
  %  dim  - dimension of intersection

  if length(rec1)>1 || length(rec2)>1
    error('Can not run isNeighbor on arrays')
  end

  ret = false; dim = [];
  isect = intersect(rec1, rec2);
  if length(isect.getFlatDims) == 1
    dim = isect.getFlatDims;
    ret = true;
  end
end