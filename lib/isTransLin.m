function result = isTransLin(rec1, rec2, vField)

% Function to check if there exists a transition from some point in rec1
% to some point in rec2 under the _linear_ vector field vField.
% The computation is done by evaluating the normal flow at corner points.
% Return true if there exists a flow and false if one can 
% find a certificate that guarantees the non-existence.

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
    h1 = zeros(1,rec1.dim);
    if rec1.xmax(flatdim)==rec2.xmin(flatdim)
        h1(flatdim)=1;
    else
        h1(flatdim)=-1;
    end 

    result = false;
    n = length(irec.getFullDims); % number of remaining dimensions
    for i = 1:2^n % iterate over vertices
        vert = irec.getVertexI(i);
        if h1*(vField.A*vert'+vField.K) > 0
            result = true;  % there is a flow
            return;
        end
    end
end