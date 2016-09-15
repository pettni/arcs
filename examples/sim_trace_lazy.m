function [xvec avec] = sim_trace_lazy(part, winning, cinv, K, act_set, trans)
    maxdisctrans = 100;
    dt = 1;
    if length(winning) > length(cinv)
        winningnotcinv = setdiff(winning, cinv);
    else
        winningnotcinv = winning;
    end

    dim = length(act_set{1}.A);
    s0 = winningnotcinv(randi(length(winningnotcinv))); % need to start in winning
    rec0 = part(s0);
    x0 = rec0.getMidpoint()';

    x = x0;
    xvec = zeros(dim,0);
    avec = zeros(1,0);

    for disc_trans=1:maxdisctrans
        s = find_cell(part, x);
        if disc_trans>1;
            next_with_act = find(trans(s,:,act)==1);
            if ~isempty(setdiff(next_with_act,cinv))
                act = find_K(s, winning, K);
            end
        else
            act = find_K(s, winning, K);
        end
        rec = part(s);

        pos0 = max(1, size(xvec,2));
        while isInside(rec, x)
            dx = act_set{act}.A*x+act_set{act}.K;
            x = x+dx*dt;
            xvec(:, end+1) = x;
            avec(:, end+1) = act;
        end

    end
end

function num = find_cell(part, x0)
	num = -1;
	for i = 1:length(part)
		if part.cell_list(i).isInside(x0);
			num = i;
			return;
		end
	end
end

function k = find_K(s, winning, K)
	ind = find(winning==s);
	k = K(ind);
end