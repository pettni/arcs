clear all;
domain = Rec([20 28; 20 28; 20 28]);
goal_set = Rec([21 27; 22 25; 22 25], {'SET'});
% unsafe_set = Rec([27.7 28; 27 28; 27.2 28], {'UNSAFE'});

maxiter = 100;
split_goal = false;   % to avoid zeno

load('radiant_data/a1.mat')
load('radiant_data/a2.mat')
load('radiant_data/b1.mat')
load('radiant_data/b2.mat')

b1(3) = b1(3)/10; 
b2(3) = b2(3)/10;

tic

% Build initial partition
part = Partition(domain);
part.add_area(goal_set);

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
part.create_ts();
part.add_mode({a1, b1});
part.add_mode({a2, b2});

%%%%%%%%%%%%%%%%%%%%%%%%
% SYNTHESIS-REFINEMENT %
%%%%%%%%%%%%%%%%%%%%%%%%

Win = [];
V1 = [];
iter = 0;
while true

  if isempty(Win)
    vol = 0;
  else
    vol = sum(volume(part(Win)))/volume(domain);
  end

  N = length(part);
  
  time = toc;
  disp(['iteration ', num2str(iter), ', time ', num2str(time), ', states ', num2str(N), ', winning set volume ', num2str(vol)])

  % Solve [] A && <>[] B &&_i []<> C_list{i}
  A = 1:N;
  B = part.get_cells_with_ap({'SET'});
  C_list = {1:N};

  % Winning set
  [Win, Vlist_t] = part.ts.win_primal(A, B, C_list, 'exists', Win);
  
  if length(Vlist_t) > 0
    V1 = Vlist_t{1};
  end

  % Candidate set
  C = union(setdiff(part.ts.pre(Win, 1:2, 'exists', 'exists'), Win), ...
        setdiff(B, V1));

  if isempty(C) || iter == maxiter
    % we are done
    break
  end

  % Split largest cell in candidate set
  [~, C_index] = max(volume(part(C)));
  part.split_cell(C(C_index));

  iter = iter + 1;
end

% Get control strategy
[~, Vlist, Klist] = part.ts.win_primal(A, B, C_list, 'exists');

toc
