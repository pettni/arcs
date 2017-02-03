% clear all;
yalmip clear;

global ops
ops = sdpsettings('solver', 'mosek', 'cachesolvers', 1, 'verbose', 0);

domain = Rec([-2 -1.5; 2 3]);
goal_set = Rec([-1 0.5; -0.2 1.8], {'goal'});
unsafe_set = Rec([-2 -1.5; -.5 -1], {'unsafe'});

maxiter = 80;
show_plot = 1;
% use_pgs = 0;
% disturbance = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Build initial partition
part = Partition(domain);
part.add_area(goal_set);
part.add_area(unsafe_set)
part.check();   % sanity check

sdpvar x1 x2 d1

if ~disturbance
  fx1 = [-x2-0.5*x1^3-1.5*x1;...
      -x2^2+2+x1];
  fx2 = [-x2-0.5*x1^3-1.5*x1;...
      -x2+x1];
  fx3 = [-x2-0.5*x1^3-1.5*x1+2;...
      10+x1];
  fx4 = [-x2-0.5*x1^3-1.5*x1-1.5;...
      -10+x1];

  dyn_list = {{fx1, [x1; x2]}, ...
            {fx2, [x1; x2]}, ...
            {fx3, [x1; x2]}, ...
            {fx4, [x1; x2]}};
else
  fx1 = [-x2-0.5*x1^3-1.5*x1;...
      -x2^2+2+x1+d1];
  fx2 = [-x2-0.5*x1^3-1.5*x1;...
      -x2+x1+d1];
  fx3 = [-x2-0.5*x1^3-1.5*x1+2;...
      10+x1+d1];
  fx4 = [-x2-0.5*x1^3-1.5*x1-1.5;...
      -10+x1+d1];

  d_rec = Rec([-0.5, 0.5]);
  dyn_list = {{fx1, [x1; x2], [d1], d_rec}, ...
              {fx2, [x1; x2], [d1], d_rec}, ...
              {fx3, [x1; x2], [d1], d_rec}, ...
              {fx4, [x1; x2], [d1], d_rec}};
end

% Create transition system
part.create_ts();

% Add modes
part.add_mode(dyn_list{1})
part.add_mode(dyn_list{2});
part.add_mode(dyn_list{3});
part.add_mode(dyn_list{4});

% Disable progress groups
if ~use_pgs
  part.ts.b_disable_pg = true;
else
  % Search for transient areas
  part.search_trans_reg(3);
end

Win = [];
iter = 0;
tic
while true

  U = part.get_cells_with_ap({'unsafe'});
  A = setdiff(1:length(part), U);
  B = part.get_cells_with_ap({'goal'});

  [Win, Cwin] = part.ts.win_primal(A, B, [], 'exists', 'forall', Win);
  part.add_aps(Win, {'win'});

  Cwin = setdiff(Cwin, union(Win, U));

  time = toc;
  disp(['iteration ', num2str(iter), ', time ', num2str(time), ', states ', num2str(length(part)), ', winning set volume ', num2str(sum(volume(part(Win)))/volume(part.domain))])

  if isempty(Cwin) || iter == maxiter
    % we are done
    break
  end

  if show_plot
    clf; hold on
    part.plot_vf();
    plot(part(Cwin), [0 0 1], 0.8, 0);
    drawnow;
  end

  % Split largest cell in candidate set
  [~, C_index] = max(volume(part.cell_list(Cwin)));
  part.split_cell(Cwin(C_index));

  iter = iter + 1;
end

if show_plot
  clf; hold on
  part.plot();
  drawnow;
end

[~, ~, cont] = part.ts.win_primal(A, B, [], 'exists', 'forall');
