function ret = isTransOutLin(rec, dom, fx)
	% Assuming that rec \subset dom, see if there is a transition from
	% rec to the outside of dom.
	dim = rec.dim;
    rest = mldivide( Rec([-Inf*ones(1,dim); Inf*ones(1,dim)]), dom);
    ret = false;
    for part=rest
    	if intersects(rec, part)
	        if isTransLin(rec, part, fx)
	            ret = true;
	            return;
	        end
	    end
    end
end