data = load('run_time_data.txt');
% data(1, :) log encoding run times
% data(2, :) split encoding run times
% data(3, :) log encoding run times with reordering
% data(4, :) split encoding run times with reordering

semilogy(data(1, :), '-o');
hold on
semilogy(data(2, :), '-o');
semilogy(data(3, :), '-o');
semilogy(data(4, :), '-o');
hold off
xlabel('Test id');
ylabel('win primal CPU time (s)');
legend('Log', 'Split', 'Reordered log', 'Reordered split', ...
      'Location', 'northwest');
grid on;