function ret = pre(states, trans_set, act)
	% PRE: computes the 1-step Pre set of a set of states under the action(s) act
	% if no action is specified, the Pre set of all action is computed
	if nargin<3
		act = 1:size(trans_set,3);
	end
	ret = [];
	for a = act
		ret = union(ret, find(sum(trans_set(:,states,a),2)>=1));
	end
	ret = reshape(ret, 1, length(ret));
end