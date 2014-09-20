function [trans, trans_out] = update_transitions(fx, part, chInd, oldtrans, oldtrans_out)
	% Return new transition matrices
	% Assumes a new state has been added by splitting a state
	% The new positions are at chInd and at the end of the part vector
		
	if length(fx)>1
		trans = zeros(length(part),length(part),length(fx));
		trans_out = zeros(length(part),1,length(fx));
		for i=1:length(fx);
			[trans(:,:,i) trans_out(:,:,i)] = ...
					update_transitions(fx{i}, part, chInd, oldtrans(:,:,i), oldtrans_out(:,:,i));
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
		if isTransLin(part(chInd), part(adj), fx)
			trans(chInd, adj) = 1;
		end
		if isTransLin(part(adj), part(chInd), fx)
			trans(adj, chInd) = 1;
		end
	end

	% Adjacencies of second cell
	for adj=adj2
		if isTransLin(part(N), part(adj), fx)
			trans(N, adj) = 1;
		end
		if isTransLin(part(adj), part(N), fx)
			trans(adj, N) = 1;
		end
	end

	% Transient - if previous was transient, so are these, so only need to 
	% check non-transient case
	if oldtrans(chInd, chInd)
		trans(chInd, chInd) = not(isTransientLin(part(chInd),fx));
		trans(N, N) = not(isTransientLin(part(N),fx));
	end

	% Transitions to outside
	if oldtrans_out(chInd)
		if isTransOutLin(part(chInd), part.domain, fx)
            trans_out(chInd) = 1;
        end
        if isTransOutLin(part(N), part.domain, fx)
            trans_out(N) = 1;
        end
	end
end