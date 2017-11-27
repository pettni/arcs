matlab_data = load('matlab_test_data.txt');
% matlab_data(1, :) number of transitions in each tests
% matlab_data(2, :) runtimes of matlab implem.

bdd_data = load('run_time_data.txt');
% bdd_data(1, :) log encoding run times
% bdd_data(2, :) split encoding run times
% bdd_data(3, :) log encoding run times with reordering
% bdd_data(4, :) split encoding run times with reordering

[trans, I] = sort(matlab_data(1,:));
semilogy(trans, bdd_data(1, I), '-o');
hold on
semilogy(trans, bdd_data(2, I), '-o');
semilogy(trans, bdd_data(3, I), '-o');
semilogy(trans, bdd_data(4, I), '-o');
semilogy(trans, matlab_data(2, I), '-o');
hold off
xlabel('Transition count');
ylabel('win primal CPU time (s)');
legend('Log', 'Split', 'Reordered log', 'Reordered split', ...
    'List representation', 'Location', 'NorthEastOutside');
grid on;