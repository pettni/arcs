function [V, Vlist, Klist] = win_primal(ts, A, B, C_list, quant1, V)
    % Compute winning set of
    %  []A && <>[]B &&_i []<>C_i
    % under (quant1, forall)-controllability
    % with the initial condition V ("warm start")
    
    if nargin<6
        V = uint32([]);
    end

    V = uint32(V);

    ts.create_fast();

    Vlist = {};
    Klist = {};

    while true
        Z = ts.pre(V, 1:ts.n_a, quant1, 'forall');
        for i=1:length(ts.pg_U)
            Z = union(Z, ...
                      ts.pginv(ts.pg_U{i}, ts.pg_G{i}, V, A, quant1));
        end
        [Vt, Kt] = ts.win_intermediate(A, B, Z, C_list, quant1);

        if length(Vt) == length(V)
            break
        end

        Vlist{end+1} = Vt;
        Klist{end+1} = Kt;

        V = Vt;
    end
end
