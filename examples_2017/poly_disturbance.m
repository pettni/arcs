clear all;
yalmip clear;

global ops
ops = sdpsettings('solver', 'mosek', 'cachesolvers', 1, 'verbose', 0);

domain = Rec([-2 -1.5; 2 3]);
goal_set1 = Rec([0.5 1; 1.5 2], {'goal1'});
goal_set2 = Rec([-1.5 -1.5; -1 -1], {'goal2'});
unsafe_set = Rec([-0.5 0; 0.5 1], {'unsafe'});

maxiter = 80;
show_plot = 1;
use_pgs = 1;

% Disturbance bound
d_rec = Rec([-5; 5]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sdpvar x1 x2 d1

%Moore Greitzer
fx1 = [-x2-0.5*x1^3-1.5*x1;...
    -x2^2+2+x1+d1];
fx2 = [-x2-0.5*x1^3-1.5*x1;...
    -x2+x1+d1];
fx3 = [-x2-0.5*x1^3-1.5*x1+2;...
    10+x1+d1];
fx4 = [-x2-0.5*x1^3-1.5*x1-1.5;...
    -10+x1+d1];

% Build initial partition
part = Partition(domain);
part.add_area(goal_set1);
part.add_area(goal_set2);
part.add_area(unsafe_set)
part.check();   % sanity check

% Create transition system
part.create_ts();

% Disable progress groups
if ~use_pgs
  part.ts.b_disable_pg = true;
end

% Add modes
part.add_mode({fx1, [x1; x2], [d1], d_rec});
part.add_mode({fx2, [x1; x2], [d1], d_rec});
part.add_mode({fx3, [x1; x2], [d1], d_rec});
part.add_mode({fx4, [x1; x2], [d1], d_rec});

Win = [];
iter = 0;
tic

while true

  U = part.get_cells_with_ap({'unsafe'});
  A = setdiff(1:length(part), U);
  B = part.get_cells_with_ap({'goal'});

  [Win, Cwin] = part.ts.win_primal(A, B, [], 'exists', Win);
  part.add_aps(Win, {'win'});

  Cwin = setdiff(Cwin, union(Win, U));

  time = toc;
  disp(['iteration ', num2str(iter), ', time ', num2str(time), ', states ', num2str(length(part)), ', winning set volume ', num2str(sum(volume(part(Win)))/volume(part.domain))])

  if show_plot
    part.plot();
    drawnow;
  end

  if isempty(Cwin) || iter == maxiter
    % we are done
    break
  end

  % Split largest cell in candidate set
  [~, C_index] = max(volume(part(Cwin)));
  part.split_cell(Cwin(C_index));

  iter = iter + 1;
end
