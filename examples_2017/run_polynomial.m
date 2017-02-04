for use_pgs = [0, 1]
	for disturbance = [0, 1]
		polynomial
    dstring = datestr(now);
    save(strcat('save', dstring, '.mat'), 'part', 'Win', 'cont')
	end
end