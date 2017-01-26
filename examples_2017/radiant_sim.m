% simulate radiant example using Euler forward

numsim = 50;
maxT = 48 * 3600;  % seconds
dt = 15;            % seconds

xvec_list = {};
avec_list = {};

act_set = {{a1, k1}, {a2, k2}};

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
  while isInside(boringset, x)
    x = win_min' + rand(3,1) .* (win_max' - win_min');
  end
  s = part.find_cell(x);
  t = 0;

  for i=length(Klist):-1:1
    if ismember(s, Vlist{i})
      k_counter = i;
    end
  end

  acts = Klist{k_counter}{1}(s);
  act = acts(1);

  xvec = zeros(3,ceil(maxT/dt));
  avec = zeros(1,ceil(maxT/dt));

  for i=1:ceil(maxT/dt)

    if ~isInside(part(s), x)
      % Discrete state changed
      s = part.find_cell(x);

      % Did k_counter fall?
      if ismember(s, Vlist{max(k_counter-1, 1)})
        k_counter = max(k_counter-1, 1);
      end

      % Lazy strategy
      acts = Klist{k_counter}{1}(s);
      if ~ismember(act, acts)
        act = acts(randi(length(acts)));
      end
    end

    % Sanity checks
    assert(isInside(part(s), x))
    assert(ismember(s, Vlist{k_counter}))
    x = x + (act_set{act}{1} * x + act_set{act}{2}) * dt;
    t = t + dt;
    xvec(:, i) = x;
    avec(:, i) = act;
  end

  xvec_list{end+1} = xvec;
  avec_list{end+1} = avec;

end

all_in=1;
for i=1:numsim
  if ~isInside(goal_set, xvec_list{i}(:, end))
    disp(['number ', num2str(i), ' outside'])
    all_in = 0
  end
end

disp(['all inside: ', num2str(all_in)])

save('radiant_data/radiant_sol.mat', 'part', 'Win', 'xvec_list', 'avec_list', 'goal_set', 'dt')