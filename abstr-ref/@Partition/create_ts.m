function create_ts(part, sys_set, enc_set)
	% create_ts: create a transition system representing abstraction
	if nargin < 2
    sys_set = TransSyst.sparse_set;
  end
  if nargin < 3
    enc_set = TransSyst.split_enc;
  end
  part.ts = TransSyst(length(part)+1, 0, sys_set, enc_set);
