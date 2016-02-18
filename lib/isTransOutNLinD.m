function ret = isTransOutNLinD(rec, dom, fx, vars, dBound)
	% Assuming that rec \subset dom, see if there is a transition from
	% rec to the outside of dom.
	dim = rec.dim;
    rest = mldivide( Rec([-Inf*ones(1,dim); Inf*ones(1,dim)]), dom);
    ret = false;
    for part=rest
    	if intersects(rec, part)
	        if isTransNLinD(rec, part, fx,vars, dBound)
	            ret = true;
	            return;
	        end
	    end
    end
end