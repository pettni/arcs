clear all;
global ops
global opt_settings

ops = sdpsettings('solver', 'mosek', 'cachesolvers', 1, 'verbose', 0);
opt_settings.mode = 'sdsos';
opt_settings.max_deg = 4;

system_setting = TransSyst.bdd_set;
encoding_setting = BDDSystem.log_enc;

% max # of synthesis-refinement steps
maxiter = 10000;

% Target set
goal_set = Rec([21 27; 22 25; 22 25], {'SET'});

% split final invariant set further to avoid zeno
split_inv = true;

% Disturbance: unit is W/m^2 --- heat added per unit floor area
dmax = 0.;

% Which model to use ('cdc2014' or 'realistic')
model = 'realistic';

% Progress group search depth
pg_depth = 2;

%%%%%%%%%%%%%%%%%%%%%%%
% Initial abstraction %
%%%%%%%%%%%%%%%%%%%%%%%

% Load model
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
  part.abstract({fx1, fx2}, [xvar; dvar], d_rec, system_setting, encoding_setting);
else
  fx1 = a1 * xvar + k1;
  fx2 = a2 * xvar + k2;
  part.abstract({fx1, fx2}, [xvar], [], system_setting, encoding_setting);
end

% Search for transient regions
part.search_trans_reg(pg_depth);

%%%%%%%%%%%%%%%%%%%%%%%%
% SYNTHESIS-REFINEMENT %
%%%%%%%%%%%%%%%%%%%%%%%%

Win = [];
iter = 0;
max_time = 3600;
data = zeros(maxiter, 6);
elapsed_time = 0;
part.ts.bdd_sys.dyn_reordering(false);
split_time = 0;
prime_time = 0;
%load('part.mat');
while true
  
  time = toc;
  data(iter+1, :) = [iter, time, prime_time, split_time, length(part), sum(volume(part.cell_list(Win)))/volume(part.domain)];
  disp(['iteration ', num2str(iter), ', time ', num2str(time), ', states ', num2str(length(part)), ' winning set volume: ', num2str(sum(volume(part.cell_list(Win)))/volume(part.domain))]);
  disp(['iter: ', num2str(iter), ', prime_time: ', num2str(prime_time), ' | split_time: ', num2str(split_time)]);
  % Solve <>[] 'SET'
  B = part.get_cells_with_ap({'SET'});
  
%   if mod(iter,100) == 0
%     disp('come here');
%     sys_nodes = part.ts.bdd_sys.count_system_nodes();
%     nodes = part.ts.bdd_sys.count_nodes();
%     fprintf('nodes: %d, system nodes: %d\n', nodes, sys_nodes);
%     %part.ts.bdd_sys.reorder(1000);
%     sys_nodes = part.ts.bdd_sys.count_system_nodes();
%     nodes = part.ts.bdd_sys.count_nodes();
%     fprintf('new nodes: %d, new system nodes: %d\n', nodes, sys_nodes);
%     disp('came here');
%   end
  
  prime_time = toc;
  [Win, Cwin] = part.ts.win_primal([], B, [], 'exists', 'forall', Win);
  prime_time = toc - prime_time;
  % No need to split inside winning set
  Cwin = setdiff(Cwin, Win);

  if isempty(Cwin) || iter == maxiter || time >= max_time
    break;
  end

%   Split largest cell in candidate set
  [~, C_index] = max(volume(part.cell_list(Cwin)));
  split_time = toc;
  part.split_cell(Cwin(C_index));
  split_time = toc - split_time;
  
  iter = iter + 1;
end

save('end_to_end_split.mat', 'data');
% Split final set to eliminate Zeno
% if split_inv
%   inv_set = part.ts.win_primal(part.get_cells_with_ap({'SET'}), ...
%                                [], [], 'forall');
%   for i=1:6^3
%     [~, C_index] = max(volume(part(inv_set)));
%     [ind1, ind2] = part.split_cell(inv_set(C_index));
%     inv_set = union(inv_set, ind2);
%   end
% end
% 
% Get control strategy
% [~, ~, cont] = part.ts.win_primal([], part.get_cells_with_ap({'SET'}), [], 'exists', 'forall');

toc
