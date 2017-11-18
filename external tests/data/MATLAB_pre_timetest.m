transition_file = 'BDD_statetest5_log.txt';
test_file = '2000_t_2_a_500_s.txt';

state_num = 0;
transition_num = 0;

in = [];
out = [];
action = [];

fileID = fopen(transition_file, 'r');
input = fscanf(fileID, '%*c %d %d\n', [1, 2]);
state_num = input(1);
input = fscanf(fileID, '%d %s\n', [state_num+1, 2]);
input = fscanf(fileID, '%*c %d %d\n', [1, 2]);
action_num = input(1);
transition_num = input(2);

in = zeros(transition_num, 1);
out = zeros(transition_num, 1);
action = zeros(transition_num, 1);

for i = 1:transition_num
    in(i, 1) = fscanf(fileID, '%d', [1,1]);
    a_enc = fscanf(fileID, '%s', [1, 1]);
    action(i) = bin2dec(a_enc) + 1;
    out(i, 1) = fscanf(fileID, '%d', [1,1]);
end

fclose(fileID);

ts = TransSyst(state_num+1, action_num);
for i = 1:transition_num
    ts.add_transition(in(i), out(i), action(i));
end
ts.create_fast();

fileID = fopen(test_file, 'r');
input = fscanf(fileID, '%d', [3, 1]);
test_num = input(3);

test_states = {};
test_actions = {};
test_ids = {};

for i = 1:test_num
    input = fscanf(fileID, '%*s %d %*s %d', [1, 2]);
    test_state_num = input(1);
    test_action_num = input(2);
    test_ids{i} = fscanf(fileID, '%s', [1, 2]);
    test_states{i} = fscanf(fileID, '%d', [1, test_state_num]);
    test_actions{i} = fscanf(fileID, '%d', [1, test_action_num]);
end

fclose(fileID);

pre_ee_times = [];
pre_ea_times = [];
pre_ae_times = [];
pre_aa_times = [];

for i = 1:test_num
    ids = test_ids{i} == 'e';
    tic;
    answer = ts.pre(test_states{i}, test_actions{i}, ids(1), ids(2));
    time = toc;
    if (ids(1) && ids(2))
        pre_ee_times(end+1) = time;
    elseif (ids(1) && ~ids(2))
        pre_ea_times(end+1) = time;
    elseif (~ids(1) && ids(2))
        pre_ae_times(end+1) = time;
    elseif (~ids(1) && ~ids(2))
        pre_aa_times(end+1) = time;
    end
end

ee_avtime = sum(pre_ee_times);
ea_avtime = sum(pre_ea_times);
ae_avtime = sum(pre_ae_times);
aa_avtime = sum(pre_aa_times);
disp(['Times measured: ee = ', num2str(ee_avtime), ' ea = ', num2str(ea_avtime), ...
    ' ae = ', num2str(ae_avtime), ' aa = ', num2str(aa_avtime)]);
disp(['Total time: ', num2str(sum([pre_ee_times, pre_ea_times, pre_ae_times, pre_aa_times]))]);



