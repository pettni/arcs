function result = isTransNLin(rec1, rec2, act, deg)

% Function to check if there exists a transition from some point in rec1
% to some point in rec2 under the _linear_ vector field fx.
% The computation is done by evaluating the normal flow at corner points.
% Return true if there exists a flow and false if one can 
% find a certificate that guarantees the non-existence.

    fx = act{1};
    vars = act{2};

    irec = intersect(rec1,rec2);
    fx = fx{1};
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
    h1 = zeros(1,rec1.dim);
    if rec1.xmax(flatdim)==rec2.xmin(flatdim)
        h1(flatdim)=1;
    else
        h1(flatdim)=-1;
    end 
    % h1 normal vector

    result = true;
    n = length(irec.getFullDims); % number of remaining dimensions
    poly = -h1*(fx);
    F = [];
    for i=1:length(vars)
        F = [F; irec.xmin(i)<=vars(i)<=irec.xmax(i)];
    end
    global ops;
    diagnosis = solvemoment(F, poly, ops, deg);
    if (diagnosis.problem==0 || diagnosis.problem==4) & value(poly)>0
        result = false;  % there is no flow
    end
end