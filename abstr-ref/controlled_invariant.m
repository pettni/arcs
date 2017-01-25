function [C,K] = controlled_invariant(X, trans_set)
% PRE_EXISTS_FORALL: Returns the maximal controlled-invariant set C contained in X, along
% with a contoller K that makes C invariant when it is used.
%
% SYNTAX
% ------
%
%	[C, K] = controlled_invariant(X, trans_set)
% 
% INPUT
% -----
%	
%	X 			initial set of states
%	trans_set	tensor of transitions of size (n x n x m), where m is the number of actions
%
% OUTPUT
% -----
%	
%	C 	maximal controlled-invariant set contained in X
%	K	controller that makes C controlled-invariant
%

C = X;
while 1 % termination guaranteed - finite
	[Ct, K] = pre_exists_forall(C, trans_set);
	if isempty(setdiff(C, Ct))
		% remove state/controller combinations which are not in C		
		K(~ismember(Ct, C)) = [];
		Ct(~ismember(Ct, C)) = [];
		C = Ct;
		break;
	else
		C = intersect(Ct, C);
	end
end