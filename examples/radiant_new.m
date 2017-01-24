clear all;
addpath('../lib')
global ops
ops = sdpsettings('solver','sedumi','cachesolvers',1,'verbose',0);
tic

domain = Rec([20 28; 20 28; 20 28]);
goal_set = Rec([21 27; 22 25; 22 25], {'SET'});
% unsafe_set = Rec([27.7 28; 27 28; 27.2 28], {'UNSAFE'});

maxiter = 100;
split_goal = true;   % to avoid zeno

load('radiant_data/a1.mat')
load('radiant_data/a2.mat')
load('radiant_data/b1.mat')
load('radiant_data/b2.mat')

b1(3) = b1(3)/10; 
b2(3) = b2(3)/10;

fx1.A = a1;
fx1.K = b1;

fx2.A = a2;
fx2.K = b2;
act_set={fx1,fx2};

tic

% Build initial partition
part = Partition(domain);
part.add_area(goal_set);
% part.add_area(unsafe_set)

% Split goal set to reduce Zeno
if split_goal
	for i=1:128
		goal = part.get_cells_with_ap('SET');
		[~, C_index] = max(volume(part(goal)));
		split_index = goal(C_index);
		part.split_cell(split_index);
	end
end

part.check();   % sanity check

% Build transition system
N = length(part);
ts = TransSyst(N+1, 2);  % state N+1 is outside, required for AFTS to be welldef

% Add transitions
for act = 1:2
	for i=1:N
		% Neighbor transitions
	    [adj, ~] = part.get_neighbors(i);
	    for j=adj
	        if isTransLin(part(i), part(j), act_set{act})
	            ts.add_transition(i, j, act);
	        end
	    end

	    % Out-of-domain
	    if isTransOutLin(part(i), part.domain, act_set{act})
	        ts.add_transition(i, N+1, act);
	    end

	    % Self transitions
	    if ~isTransientLin(part(i),act_set{act})
	    	ts.add_transition(i, i, act);
	    end
	end
end

% Add progress groups
calG = cell(1,length(act_set));
for act_ind = 1:2
    calG{act_ind} = {};

    % Whole domain progress groups
    if isTransientLin(part.domain,act_set{act_ind})
    	ts.add_progress_group([act_ind], 1:N)
    else
    	% Individual state progress group
        for st_ind = 1:N
            if isTransientLin(part(st_ind), act_set{act_ind})
            	ts.add_progress_group([act_ind], [st_ind])
            end
        end
    end
end

%%%%%%%%%%%%%
% SYNTHESIS %
%%%%%%%%%%%%%

Win = [];
V1 = [];
iter = 0;
while true 	

	if isempty(Win)
		vol = 0;
	else
		vol = sum(volume(part(Win)))/volume(domain);
	end
	
	time = toc;
	disp(['iteration ', num2str(iter), ', time ', num2str(time), ', states ', num2str(N), ', winning set volume ', num2str(vol)])

	N = length(part);

	% Want [] A && <>[] B
	A = 1:N;
	% A = setdiff(1:N, part.get_cells_with_ap({'UNSAFE'}));
	B = part.get_cells_with_ap({'SET'});
	C_list = {1:N};

	% Winning set
	[Win, Vlist_t] = ts.win_primal(A, B, C_list, 'exists', Win);
	
	if length(Vlist_t) > 0
		V1 = Vlist_t{1};
	end

	% Candidate set
	C = union(setdiff(ts.pre(Win, 1:2, 'exists', 'exists'), Win), ...
	          setdiff(B, V1));
	if isempty(C) || iter == maxiter
		break
	end

	% Split largest cell in candidate set
	[~, C_index] = max(volume(part(C)));
	split_index = C(C_index);
	part.split_cell(split_index);

	% Update TS
	ts.split_state(split_index);

	% Update transitions
	for act=1:2
		for i=[split_index, N+1]
		    [adj, dim] = part.get_neighbors(i);
		    for j=adj
		        if isTransLin(part(i), part(j), act_set{act})
		            ts.add_transition(i, j, act);
		        end
		        if isTransLin(part(j), part(i), act_set{act})
		            ts.add_transition(j, i, act);
		        end
		    end

		    % Out-of-domain
		    if isTransOutLin(part(i), part.domain, act_set{act})
		        ts.add_transition(i, N+2, act);
		    end

		    % Self transitions
		    if ~isTransientLin(part(i),act_set{act})
		    	ts.add_transition(i, i, act);
		    end
		end
	end

	iter = iter + 1;
end

% Get control strategy
[~, Vlist, Klist] = ts.win_primal(A, B, C_list, 'exists');

toc