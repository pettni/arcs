function result = is_trans_nlin(rec1, rec2, dyn, drec, deg)
  % is_trans_nlin(rec1, rec2, dyn): Return false if a certificate that 
  % guarantees the non-existence of a flow from rec1 to rec2 under 
  % nonlinear dynamics dyn is found, true otherwise
  % 
  % Inputs 
  %   - rec1, rec2: sets
  %   - dyn = {fx, vars, drec}: dynamics \dot x = fx(x,d), d \in drec
  %   - drec: disturbance Rec
  %   - deg: relaxation order for moment hierarchy

  global ops;
  result = true;

  fx = dyn{1};
  vars = dyn{2};

  if isempty(drec)
    drec = Rec(zeros(0,2));
  end

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
  
  % Combined bounds
  crec = irec * drec;
	F = [crec.xmin' <= vars, vars <= crec.xmax'];

	diagnosis = solvemoment(F, obj, ops, deg);
	if (diagnosis.problem==0 || diagnosis.problem==4) & value(obj)>0
		result = false;  % there is no flow
	end
end