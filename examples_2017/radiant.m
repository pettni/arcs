clear all;
maxiter = 650;
split_inv = true;   % to avoid zeno

goal_set = Rec([21 27; 22 25; 22 25], {'SET'});

[a1 k1 e1 a2 k2 e2] = radiant_dyn();

% Disturbance: unit is W/m^2 --- heat added per unit floor area
dmax = 3;
d_rec = Rec([-dmax -dmax; dmax dmax]);

act_set = {{a1, k1, e1, d_rec}, {a2, k2, e2, d_rec}};

tic

% Build initial partition
part = Partition(Rec([20 28; 20 28; 20 28]));
part.add_area(goal_set);
part.check();   % sanity check

% Build transition system
part.create_ts();
part.add_mode(act_set{1});
part.add_mode(act_set{2});

%%%%%%%%%%%%%%%%%%%%%%%%
% SYNTHESIS-REFINEMENT %
%%%%%%%%%%%%%%%%%%%%%%%%

Win = [];
iter = 0;
while true

  if isempty(Win)
    vol = 0;
  else
    vol = sum(volume(part(Win)))/volume(part.domain);
  end

  time = toc;
  disp(['iteration ', num2str(iter), ', time ', num2str(time), ', states ', num2str(length(part)), ', winning set volume ', num2str(vol)])

  % Solve <>[] 'SET'
  [Win, Cwin] = part.ts.win_primal([], part.get_cells_with_ap({'SET'}), ...
                                   [], 'exists', Win);

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

% Split final set to eliminate Zeno
if split_inv
  inv_set = part.ts.win_primal(part.get_cells_with_ap({'SET'}), ...
                               [], [], 'exists');
  for i=1:6^3
    [~, C_index] = max(volume(part(inv_set)));
    [ind1, ind2] = part.split_cell(inv_set(C_index));
    inv_set = union(inv_set, ind2);
  end
end

% Get control strategy
[Win, ~, Vlist, Klist] = ...
    part.ts.win_primal([], part.get_cells_with_ap({'SET'}), [], 'exists');

toc
