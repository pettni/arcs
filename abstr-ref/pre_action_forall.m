function pre_set = pre_action_forall(X, trans_set, act)
	% Computes the set pre_set, such that for all s \in pre_set,
	% s -> X under action act.
	N_states = size(trans_set,2);
	pre_cand = pre(X, trans_set, act);
	bad_states = setdiff(post(pre_cand, trans_set, act), X);
	pre_set = pre_cand(find(sum(trans_set(pre_cand,bad_states,act),2) == 0));
	pre_set = reshape(pre_set, 1, length(pre_set));
end