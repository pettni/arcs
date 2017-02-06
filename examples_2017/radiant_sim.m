% simulate radiant example using Euler forward

numsim = 50;
maxT = 72 * 3600;  % seconds
dt = 20;           % seconds

xvec_list = cell(1, numsim);
avec_list = cell(1, numsim);

% Assume winning set is a square
win_min = [Inf Inf Inf];
win_max = [-Inf -Inf -Inf];
for rec = part(Win)
  win_max = max(win_max, rec.xmax);
  win_min = min(win_min, rec.xmin);
end

boringset = Rec([win_min + 0.2*(win_max - win_min);
                win_max - 0.2*(win_max - win_min)]');

randseed(1);

for sim = 1:numsim
  x = win_min' + rand(3,1) .* (win_max' - win_min');
  while isInside(boringset, x) || part.find_cell(x) == -1 || ~ismember(part.find_cell(x), Win)
    x = win_min' + rand(3,1) .* (win_max' - win_min');
  end
  s = part.find_cell(x);
  t = 0;

  acts = cont(s);
  act = acts(1);

  offset1 = rand(1,1);
  offset2 = rand(1,1);

  xvec = zeros(3,ceil(maxT/dt));
  avec = zeros(1,ceil(maxT/dt));

  for i=1:ceil(maxT/dt)

    if ~isInside(part.cell_list(s), x)
      % Discrete state changed
      s = part.find_cell(x);

      % Lazy strategy
      acts = cont(s);
      if ~ismember(act, acts)
        act = acts(randi(length(acts)));
      end
    end

    % Get disturbance
    d1 = dmax * sin(2*pi*offset1 + 2*pi*t/(24*3600));
    d2 = dmax * sin(2*pi*offset2 + 2*pi*t/(24*3600));

    x = x + (part.dyn_list{act}{1} * x + ...
             part.dyn_list{act}{2} + ...
             part.dyn_list{act}{3} * [d1; d2]) * dt;
    t = t + dt;
    xvec(:, i) = x;
    avec(:, i) = act;
  end

  xvec_list{sim} = xvec;
  avec_list{sim} = avec;

end

all_in=1;
for i=1:numsim
  if ~isInside(goal_set, xvec_list{i}(:, end))
    disp(['number ', num2str(i), ' outside'])
    all_in = 0;
  end
end

disp(['all inside: ', num2str(all_in)])

% save('radiant_data/radiant_sol_paper.mat', 'part', 'Win', 'xvec_list', 'avec_list', 'goal_set', 'dt')