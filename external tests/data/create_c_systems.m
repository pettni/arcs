function create_c_systems(id)

name = 'BDD_statetest';
s = load([name, num2str(id), '.mat']);
part = s.part;
ts = part.ts;

% split encodings
file = fopen([name, num2str(id), '_split.txt'], 'w');
N = length(part.cell_list);
var_num = 0;
for i = 1:N
    cell = part.cell_list(i);
    enc_length = length(cell.encoding);
    var_num = max([var_num, enc_length]);
end
var_num = var_num + 1;

fprintf(file, 'N %d %d\n', N, var_num);
for i = 1:N
    cell = part.cell_list(i);
    encoding = sprintf('%d', [1, cell.encoding, zeros(1, var_num - length(cell.encoding)-1)]);
    fprintf(file, '%d %s\n', i, encoding);
end
encoding = sprintf('%d', zeros(1, var_num));
fprintf(file, '%d %s\n', N+1, encoding);

T = ts.n_a;
TS = part.ts;
in = TS.state1;
out = TS.state2;
a = TS.action;
a_vars = ceil(log2(double(T)));
fprintf(file, 'T %d %d\n', ts.n_a, length(in));
for i = 1:length(in)
    encoding = sprintf('%d', de2bi(a(i)-1));
    encoding = [encoding, repmat('0', 1, -(length(encoding)-a_vars))];
    fprintf(file, '%d %s %d\n', in(i), encoding, out(i));
end 

pg_group_num = length(ts.pg_G);
fprintf(file, 'Groups %d\n', pg_group_num);
for i = 1:pg_group_num
    pg_group_size_U = length(ts.pg_U{i});
    pg_group_size_G = length(ts.pg_G{i});
    fprintf(file, 'Group_%d %d %d\n', i, pg_group_size_U, pg_group_size_G);
    fprintf(file, '%d ', ts.pg_U{i});
    fprintf(file, '%d ', ts.pg_G{i});
    fprintf(file, '\n');
end
fclose(file);

% log encodings
file = fopen([name, num2str(id), '_log.txt'], 'w');
N = length(part.cell_list);
vars = ceil(log2(N+1));
fprintf(file, 'N %d %d\n', N, vars);
for i = 1:N
    cell = part.cell_list(i);
    encoding = sprintf('%d', de2bi(i));
    encoding = [encoding, repmat('0', 1, -(length(encoding)-vars))];
    fprintf(file, '%d %s\n', i, encoding);
end
encoding = sprintf('%d', de2bi(0));
encoding = [encoding, repmat('0', 1, -(length(encoding)-vars))];
fprintf(file, '%d %s\n', N+1, encoding);

T = ts.n_a;
TS = part.ts;
in = TS.state1;
out = TS.state2;
a = TS.action;
a_vars = ceil(log2(double(ts.n_a)));
fprintf(file, 'T %d %d\n', ts.n_a, length(in));
for i = 1:length(in)
    encoding = sprintf('%d', de2bi(a(i)-1));
    encoding = [encoding, repmat('0', 1, -(length(encoding) - a_vars))];
    fprintf(file, '%d %s %d\n', in(i), encoding, out(i));
end

pg_group_num = length(ts.pg_G);
fprintf(file, 'Groups %d\n', pg_group_num);
for i = 1:pg_group_num
    pg_group_size_U = length(ts.pg_U{i});
    pg_group_size_G = length(ts.pg_G{i});
    fprintf(file, 'Group_%d %d %d\n', i, pg_group_size_U, pg_group_size_G);
    fprintf(file, '%d ', ts.pg_U{i});
    fprintf(file, '%d ', ts.pg_G{i});
    fprintf(file, '\n');
end

fclose(file);

end