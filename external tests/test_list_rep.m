%% perform test
fprintf("Starting MATLAB list representation tests\n");
addpath(genpath('data/'))
addpath(genpath('../abstr-ref'))

tests = 7;
runtime_file = 'matlab_test_data.txt';
runtimes = zeros(1, tests);
transition_counts = zeros(1, tests);

data_path = 'data/';

for i = 1:tests
  test_path = sprintf("test%d/BDD_statetest.mat", i);
  test_data = load(data_path + test_path);
  part = test_data.part;
  transition_counts(i) = part.ts.num_trans();
  
  B = part.get_cells_with_ap({'SET'});
  t = cputime;
  Win = part.ts.win_primal([], B, [], 'exists', 'forall');
  runtime = cputime - t;
  runtimes(i) = runtime;
  fprintf("Matlab runtime test %d complete\n", i);
end

%% dump runtime data to file
% warning: overwrites old data
fileID = fopen(runtime_file, 'w');
fprintf(fileID, "%d ", transition_counts);
fprintf(fileID, "\n");
fprintf(fileID, "%.6f ", runtimes);
fclose(fileID);




