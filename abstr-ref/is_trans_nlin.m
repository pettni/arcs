function result = is_trans_nlin(rec1, rec2, dyn, deg)
  % is_trans_nlin(rec1, rec2, dyn): Return false if a certificate that 
  % guarantees the non-existence of a flow from rec1 to rec2 under 
  % nonlinear dynamics dyn is found, true otherwise
  % 
  % Inputs 
  %   - rec1, rec2: sets
  %   - dyn = {fx, x, (d, drec)}: dynamics \dot x = fx(x,d), d \in drec
  %   - deg: relaxation order for moment hierarchy

  global ops;
  result = true;

  irec = intersect(rec1,rec2);
  % overlapping
  if irec.isFullDim
    warning('isTransLinRec: was called with overlapping Recs - returning true');
    result = true;
    return;
  end

  % are they adjacent?
  if length(irec.getFlatDims) ~= 1
    result = false;
    return;
  end

	% find "flat" dimension
	flatdim = irec.getFlatDims;
	h1 = zeros(1, length(dyn{1}));
	if rec1.xmax(flatdim)==rec2.xmin(flatdim)
		h1(flatdim)=1;
	else
		h1(flatdim)=-1;
	end 
	% h1 normal vector

	obj = -h1*dyn{1};
	F = [irec.xmin' <= dyn{2}, dyn{2} <= irec.xmax'];
	if length(dyn) > 2
		% Disturbance
		dvar = dyn{3};
		drec = dyn{4};
		F = [F, drec.xmin' <= dvar, dvar <= drec.xmax'];
	end

	diagnosis = solvemoment(F, obj, ops, deg);
	if (diagnosis.problem==0 || diagnosis.problem==4) & value(obj)>0
		result = false;  % there is no flow
	end
end