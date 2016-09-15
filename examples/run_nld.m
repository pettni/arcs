% Arbitrary 2-d nonlinear model

% 7/12/13: adding rectangular adjacency
% 9/19/13: convert to mpt3

%clear all;
clear all; clf;
addpath('../lib')
global ops
ops = sdpsettings('solver','sedumi','cachesolvers',1,'verbose',0);
tic

% Load the problem here!

% load_linear_codim_engine; 
% load_linear_codim_one 
% load_radiant;
load_poly_disturbance;
% load_poly_lin_test;

maxiter = 80;
show_plot = 0;

% Compute transitions
[trans, trans_out] = compute_transitions_nld(act_set, part, vars, dBound);

deg = 4;
% Compute progress group
%calG = compute_progress_group_nld(act_set, part, vars, 6, dBound);
calG = {{},{},{},{}};
%%
for iter = 1:maxiter
    iter
    N_state = length(part);
    N_act = length(act_set);

    % add dummy state representing outside
    trans_set = zeros(N_state+1, N_state+1, N_act);
    trans_set(1:N_state,N_state+1,:) = trans_out;
    trans_set(1:N_state,1:N_state,:) = trans;
    
    % get cinv/goal/winning/unsafe sets and control action from partition
    unsafe = [part.get_cells_with_ap(2) N_state+1]; % add outside 'state' as unsafe
    winning = part.get_cells_with_ap(3);
    cinv = part.get_cells_with_ap(4);

    K = zeros(1,length(winning));
    for act=1:size(trans_set,3)
        actset = part.get_cells_with_ap(10+act);
        K(ismember(winning, actset)) = act;
    end

    if isempty(winning)
        % we do not have a controlled-invariant set yet, let's try to find one contained in goal
        goal = part.get_cells_with_ap(1);
        [cinv, K] = controlled_invariant(goal, trans_set);
       
        cinv0 = cinv;
        Kcinv0 = K;
        % save control action as APs
        for act=1:size(trans_set,3)
            part.add_aps(cinv(find(K==act)), 10+act);
        end

        if ~isempty(cinv)
            % we found a controlled-invariant set, this is the winning set!
            disp(['Controlled-invariant set found in iteration ', num2str(iter)])
            part.add_aps(cinv, 4); % cinv
            part.add_aps(cinv, 3); % winning
            continue;
        else
            % we didn't find a controlled-invariant set, let's split the largest goal cell
            
            

            split_list = goal;


            % TO ADD: remove sets from goal which are in goal \ pre_forall_forall goal and not transient
        end
    else
        % expand winning set
        winning_1 = [];
        while ~isempty(setdiff(winning, winning_1));
            winning_1 = winning;
            [winning K split_list] = expand_winning(winning_1, K, trans_set, unsafe, calG);
        end

        for ssss = 1:length(cinv0)
            cinv_set = cinv0(ssss);
            k1 = Kcinv0(ssss);
            k2 = K(find(winning == cinv_set));
            if ~(k1 == k2)
                cinv0
                Kcinv0
                winning
                K
                error('asd')
            end
        end
        part.add_aps(winning, 3);

        for act=1:size(trans_set,3)
            part.add_aps(winning(find(K==act)), 10+act);
        end
        
        if (isempty(split_list)) break; end % nothing more to split, we are done
    end

    % expand losing set
    unsafe_1 = [];
    while ~isempty(setdiff(unsafe, unsafe_1))
        unsafe_1 = unsafe;
        unsafe = union(unsafe_1, pre_forall_forall(unsafe_1, trans_set));
    end
    part.add_aps(setdiff(unsafe, N_state+1), 2);
%for subst = 1:10
    %split largest cell in in split_list
    [~, maxvolind] = max(volume(part(split_list)));
    chInd = split_list(maxvolind);
    part.split_cell(chInd);

    % update transitions locally around chInd
    [trans, trans_out] = update_transitionsNLD(act_set, part, chInd, trans, trans_out,vars, dBound);

    %update calG locally around chInd
    calG = update_progress_group(calG, part, chInd);
%end
    % plot new partition
    if show_plot
        plot(part)
        drawnow
    end
end

plot(part)
set(gca,'FontSize',20)
toc
