% simualte a controller using euler forward

function simulate(part, winning, cinv, K, act_set)
	numsim = 50;
	maxdisctrans = 500;
	dt = 0.001;
	winningnotcinv = setdiff(winning, cinv);

	clf;
	hold on
	plot(part.cell_list(part.get_cells_with_ap(4)), [0 1 0], 0.1, 0);
	lims = [part.domain.xmin; part.domain.xmax];
	xlim(lims(:,1))
	ylim(lims(:,2))
	zlim(lims(:,3))

	for sim = 1:numsim
		s0 = winningnotcinv(randi(length(winningnotcinv))); % need to start in winning
		rec0 = part(s0);
		x0 = rec0.getMidpoint()';

		x = x0;
		xvec = zeros(3,0);
		avec = zeros(1,0);
	
		plot3(x0(1), x0(2), x0(3), '*k')

		for disc_trans=1:maxdisctrans
			s = find_cell(part, x);
			act = find_K(s, winning, K);
			rec = part(s);

			pos0 = max(1, size(xvec,2));
			while isInside(rec, x)
				dx = act_set{act}.A*x+act_set{act}.K;
				x = x+dx;
				xvec(:, end+1) = x;
				avec(:, end+1) = act;
			end

			color = [act==1 0 act==2];
			plot3(xvec(1,pos0:end), xvec(2,pos0:end), xvec(3,pos0:end), 'color', color)
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