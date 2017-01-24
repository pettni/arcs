% simulate radiant example using Euler forward

numsim = 20;
maxT = 1.5 * 3600;  % seconds
dt = 1;             % seconds

xvec_list = {};
avec_list = {};

start_set = setdiff(Win, B);

for sim = 1:numsim
    rec0 = part(start_set(randi(length(start_set))));

    x = rec0.getMidpoint()';
    s = part.find_cell(x);
    t = 0;

    for i=length(Klist):-1:1
        if ismember(s, Vlist{i})
            k_counter = i;
        end
    end

    acts = Klist{k_counter}{1}(s);
    act = acts(1);

    xvec = zeros(3,0);
    avec = zeros(1,0);

    while (t < maxT)

        if ~isInside(part(s), x)
            % Discrete state changed
            s = part.find_cell(x);

            % Did k_counter fall?
            if ismember(s, Vlist{max(k_counter-1, 1)})
                k_counter = max(k_counter-1, 1);
            end

            acts = Klist{k_counter}{1}(s);
            if ~ismember(act, acts)
                act = acts(1);
            end
        end

        x = x + (act_set{act}.A * x + act_set{act}.K) * dt;
        t = t + dt;
        xvec(:, end+1) = x;
        avec(:, end+1) = act;
    end

    xvec_list{end+1} = xvec;
    avec_list{end+1} = avec;

end

all_in=1;
for i=1:numsim
    if ~isInside(goal_set, xvec_list{i}(:, end))
        all_in = 0
    end
end

disp(['all inside: ', num2str(all_in)])