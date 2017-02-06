clear all;
global ops
global opt_settings

ops = sdpsettings('solver', 'mosek', 'cachesolvers', 1, 'verbose', 0);
opt_settings.mode = 'sdsos';
opt_settings.max_deg = 4;

% max # of synthesis-refinement steps
maxiter = 1000;

% split final invariant set further to avoid zeno
split_inv = true;

% Target set
goal_set = Rec([21 27; 22 25; 22 25], {'SET'});

% Disturbance: unit is W/m^2 --- heat added per unit floor area
dmax = 0.;

% Which model to use ('cdc2014' or 'realistic')
model = 'cdc2014';

% Progress group search depth
pg_depth = 2;

%%%%%%%%%%%%%%%%%%%%%%%
% Initial abstraction %
%%%%%%%%%%%%%%%%%%%%%%%

% Load model
cd radiant_data
  if strcmp(model, 'realistic')
    [a1 k1 e1 a2 k2 e2] = radiant_dyn();
  elseif strcmp(model, 'cdc2014')
    load a1; load a2; load b1; load b2
    b1(3) = b1(3)/10;
    b2(3) = b2(3)/10;
    k1 = b1;
    k2 = b2;
  else
    error('invalid model')
  end
cd ..

tic

% Build initial partition
part = Partition(Rec([20 28; 20 28; 20 28]));
part.add_area(goal_set);
part.check();   % sanity check

% Create abstraction
xvar = sdpvar(3,1);
dvar = sdpvar(2,1);

if dmax
  d_rec = Rec([-dmax -dmax; dmax dmax]);
  fx1 = a1 * xvar + k1 + e1 * dvar;
  fx2 = a2 * xvar + k2 + e2 * dvar;
  part.abstract({fx1, fx2}, [xvar; dvar], d_rec);
else
  fx1 = a1 * xvar + k1;
  fx2 = a2 * xvar + k2;
  part.abstract({fx1, fx2}, [xvar]);
end

% Search for transient regions
part.search_trans_reg(pg_depth);

%%%%%%%%%%%%%%%%%%%%%%%%
% SYNTHESIS-REFINEMENT %
%%%%%%%%%%%%%%%%%%%%%%%%

Win = [];
iter = 0;
while true

  time = toc;
  disp(['iteration ', num2str(iter), ', time ', num2str(time), ', states ', num2str(length(part)), ', winning set volume ', num2str(sum(volume(part.cell_list(Win)))/volume(part.domain))])

  % Solve <>[] 'SET'
  B = part.get_cells_with_ap({'SET'});
  [Win, Cwin] = part.ts.win_primal([], B, [], 'exists', 'forall', Win);

  % No need to split inside winning set
  Cwin = setdiff(Cwin, Win);

  if isempty(Cwin) || iter == maxiter
    break
  end

  % Split largest cell in candidate set
  [~, C_index] = max(volume(part.cell_list(Cwin)));
  part.split_cell(Cwin(C_index));

  iter = iter + 1;
end


% Split final set to eliminate Zeno
if split_inv
  inv_set = part.ts.win_primal(part.get_cells_with_ap({'SET'}), ...
                               [], [], 'forall');
  for i=1:6^3
    [~, C_index] = max(volume(part(inv_set)));
    [ind1, ind2] = part.split_cell(inv_set(C_index));
    inv_set = union(inv_set, ind2);
  end
end

% Get control strategy
[~, ~, cont] = part.ts.win_primal([], part.get_cells_with_ap({'SET'}), [], 'exists', 'forall');

toc