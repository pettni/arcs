function [pre_set, K] = pre_exists_forall(X, trans_set)
	% PRE_EXISTS_FORALL: Computes the set pre_set, such that for 
	% every s \in pre_set, there exists an action a such that s -> X under a.
	%
	% SYNTAX
	% ------
	%
	%	[pre_set, K] = pre_exists_forall(X, trans_set)
	% 
	% INPUT
	% -----
	%	
	%	X 			initial set of states
	%	trans_set	tensor of transitions (n x n x m), where m is the number of actions
	%
	% OUTPUT
	% -----
	%	
	%	pre_set 	set of states such that X can is reached under a given action
	%	K			for each state pre_set(i), K(i) is the action that guarantees that X is reached
	%

	pre_set = [];
	K = [];
	for a = 1:size(trans_set, 3)
		pre_set_a = pre_action_forall(X, trans_set, a);
		new_states = setdiff(pre_set_a, pre_set);
		pre_set = [pre_set, new_states];
		K = [K a*ones(1,length(new_states))];
	end
	pre_set = reshape(pre_set, 1, length(pre_set));
end