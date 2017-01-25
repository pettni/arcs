function [trans, trans_out] = update_transitionsNLD(fx, part, chInd, oldtrans, oldtrans_out, vars, dBound)
	% Return new transition matrices
	% Assumes a new state has been added by splitting a state
	% The new positions are at chInd and at the end of the part vector
		
	if length(fx)>1
		trans = zeros(length(part),length(part),length(fx));
		trans_out = zeros(length(part),1,length(fx));
		for i=1:length(fx);
			[trans(:,:,i) trans_out(:,:,i)] = ...
					update_transitionsNLD(fx{i}, part, chInd, oldtrans(:,:,i), oldtrans_out(:,:,i), vars, dBound);
		end
		return;
	end

	N = length(part);

	trans = zeros(N,N);
	trans_out = zeros(N,1);

	old_ind = [1:chInd-1 chInd+1:N-1];

	% keep transitions between unchanged indices, let rest be 0
	trans(old_ind, old_ind) = oldtrans(old_ind, old_ind);
	trans_out(old_ind) = oldtrans_out(old_ind);

	[adj1 dim1] = part.get_neighbors(chInd);
	[adj2 dim2] = part.get_neighbors(N);

	% Adjacencies of first cell
	for adj=adj1
		if isTransNLinD(part(chInd), part(adj), fx, vars, dBound)
			trans(chInd, adj) = 1;
		end
		if isTransNLinD(part(adj), part(chInd), fx, vars, dBound)
			trans(adj, chInd) = 1;
		end
	end

	% Adjacencies of second cell
	for adj=adj2
		if isTransNLinD(part(N), part(adj), fx, vars, dBound)
			trans(N, adj) = 1;
		end
		if isTransNLinD(part(adj), part(N), fx, vars, dBound)
			trans(adj, N) = 1;
		end
    end
    deg = 4;
	% Transient - if previous was transient, so are these, so only need to 
	% check non-transient case
	if oldtrans(chInd, chInd)
		trans(chInd, chInd) = not(isTransientNLinD(part(chInd),fx, vars, deg, dBound));
		trans(N, N) = not(isTransientNLinD(part(N),fx, vars, deg, dBound));
	end

	% Transitions to outside
	if oldtrans_out(chInd)
		if isTransOutNLinD(part(chInd), part.domain, fx, vars, dBound)
            trans_out(chInd) = 1;
        end
        if isTransOutNLinD(part(N), part.domain, fx, vars, dBound)
            trans_out(N) = 1;
        end
	end
end