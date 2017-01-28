clear all;
yalmip clear;

global ops
ops = sdpsettings('solver', 'mosek', 'cachesolvers', 1, 'verbose', 0);

domain = Rec([-2 -1.5; 2 3]);
goal_set = Rec([-1 0.5; -0.2 1.8], {'goal'});
unsafe_set = Rec([-2 -1.5; -.5 -1], {'unsafe'});

maxiter = 80;
show_plot = 1;
use_pgs = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sdpvar x1 x2
vars = [x1; x2];

%Moore Greitzer
fx1 = [-x2-0.5*x1^3-1.5*x1;...
       -x2^2+2+x1];
fx2 = [-x2-0.5*x1^3-1.5*x1;...
       -x2+x1];
fx3 = [-x2-0.5*x1^3-1.5*x1+2;...
       10+x1];
fx4 = [-x2-0.5*x1^3-1.5*x1-1.5;...
       -10+x1];

% Build initial partition
part = Partition(domain);
part.add_area(goal_set);
part.add_area(unsafe_set)
part.check();   % sanity check

% Disable progress groups
if ~use_pgs
  part.ts.b_disable_pg = true;
end

dyn_list = {{fx1, vars}, {fx2, vars}, {fx3, vars}, {fx4, vars}};

% Split some first to reveal progress groups
granularity = 1;

large_delta = (domain.xmax - domain.xmin)/granularity;
small_delta = (domain.xmax - domain.xmin)/(2*granularity);

for i = 1:2*granularity-1
    splitx = Rec([-inf -inf; domain.xmin(1) + i*small_delta(1) inf]);
    splity = Rec([-inf -inf; Inf domain.xmin(2) + i*small_delta(2)]);
    part.intersect_area(splitx);
    part.intersect_area(splity);
end

% Create transition system
part.create_ts();

% Add modes
part.add_mode(dyn_list{1})
part.add_mode(dyn_list{2});
part.add_mode(dyn_list{3});
part.add_mode(dyn_list{4});

get_pg_rec(part, Rec([-2 -1.5; 2 0.5]), 1:4);
part.ts.pg_U
get_pg_rec(part, Rec([-2 0.5; 2 3]), 1:4);
part.ts.pg_U

get_pg_rec(part, Rec([-2 -1.5; 0 3]), 1:4);
part.ts.pg_U
get_pg_rec(part, Rec([0 -1.5; 2 3]), 1:4);
part.ts.pg_U

  % for id = 0:granularity
  %   for jd = 0:granularity
  %     search_xmin = domain.xmin + [id jd].*small_delta - [0.1 0.1];
  %     search_xmax = search_xmin + large_delta + [0.1 0.1];
  %     search_rec = Rec([search_xmin; search_xmax]);
  %     get_pg_rec(part, search_rec, 1:4);
  %   end
  % end

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
  [~, C_index] = max(volume(part.cell_list(Cwin)));
  part.split_cell(Cwin(C_index));

  iter = iter + 1;
end

[Win, Cwin, Vlist, Klist] = part.ts.win_primal(A, B, [], 'exists');
