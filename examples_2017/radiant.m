clear all;
domain = Rec([20 28; 20 28; 20 28]);
goal_set = Rec([21 27; 22 25; 22 25], {'SET'});
% unsafe_set = Rec([27.7 28; 27 28; 27.2 28], {'UNSAFE'});

maxiter = 150;
split_goal = true;   % to avoid zeno

[a1 k1 a2 k2] = radiant_dyn();

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
part.add_mode({a1, k1});
part.add_mode({a2, k2});

%%%%%%%%%%%%%%%%%%%%%%%%
% SYNTHESIS-REFINEMENT %
%%%%%%%%%%%%%%%%%%%%%%%%

Win = [];
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
  [Win, Cwin] = part.ts.win_primal(A, B, C_list, 'exists', Win);

  % No need to split inside winning set
  Cwin = setdiff(Cwin, Win);

  if isempty(Cwin) || iter == maxiter
    break
  end

  % Split largest cell in candidate set
  [~, C_index] = max(volume(part(Cwin)));
  [ind1, ind2] = part.split_cell(Cwin(C_index));

  % If we happened to split winning set, update it
  if ismember(ind1, Win)
    Win(end+1) = ind2;
  end

  iter = iter + 1;
end

% Get control strategy
[~, ~, Vlist, Klist] = part.ts.win_primal(A, B, C_list, 'exists');

toc
