function pre_set = pre_forall_forall(X, trans_set)
	% Computes the set of states pre_set, such that X is reached
	% from pre_set for all actions.
	N_states = size(trans_set,2);
	pre_set = pre(X, trans_set);;
	for a = 1:size(trans_set, 3)
		pre_set = intersect(pre_set, pre_action_forall(X, trans_set, a));
	end
	pre_set = reshape(pre_set, 1, length(pre_set));
end