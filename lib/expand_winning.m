function [win1 K1 win1brd] = expand_winning(win0, K0, trans_set, unsafe, calG)
	% expands the winning set in one time step

	[S1, K1] = pre_exists_forall(win0, trans_set);
	[S1, K1] = get_new_SK(S1, K1, unsafe);						% guaranteed wins and controller
	S2 = setdiff(pre(win0, trans_set), union(S1, unsafe));		% candidate wins

	% use progress group to determine which candidates that are in fact winning
	R = S1;
	K_R = K1;
	for a = 1:size(trans_set, 3)
		S3 = [];
		for G=calG{a}
			G = setdiff(G{1}, unsafe);
			for q = intersect(S2, G)	
				S4 = [q];
				for i=1:length(G) % compute stuff reachable from q, except win0
					S4n = setdiff(union(S4, post(S4, trans_set, a)), win0);
                    if length(S4n)==length(S4)
                        break;
                    else
                        S4 = S4n;
                    end
				end
				if isempty(setdiff(S4, union(G, S1)))
					% no states except those in S1 are reachable -> Win!
					S3 = union(S3, q);
				end
			end
		end
		new_sets_a = setdiff(S3, R);
		R = [R new_sets_a];
		K_R = [K_R a*ones(1,length(new_sets_a))];
	end

	[new_sets, new_K] = get_new_SK(R, K_R, win0);
	win1 = [win0 new_sets];
	K1 = [K0 new_K];
	win1brd = setdiff(S2, win1);
end

function [S0, K0] = get_new_SK(S,K, remove)
	K0 = K;
	S0 = S;
	K0(ismember(S, remove)) = [];
	S0(ismember(S, remove)) = [];
end