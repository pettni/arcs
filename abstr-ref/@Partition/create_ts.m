function create_ts(part)
	% create_ts: create a transition system from abstraction
	part.ts = TransSyst(length(part)+1, 0);
