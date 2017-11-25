data = load('run_time_data.txt');
% data(1,:) BDD sizes in nodes (before reordering)
% data(2, :) log encoding run times
% data(3, :) split encoding run times
% data(4, :) log encoding run times with reordering
% data(5, :) split encoding run times with reordering
[size_data, I] = sort(data(1,:));

test_count = 7;

semilogy(size_data, data(2, I), '-o');
hold on
semilogy(size_data, data(3, I), '-o');
semilogy(size_data, data(4, I), '-o');
semilogy(size_data, data(5, I), '-o');
hold off
xlabel('original BDD size (in nodes)');
ylabel('win primal CPU time (s)');
legend('Log', 'Split', 'Reordered log', 'Reordered split', ...
      'Location', 'northwest');
grid on;