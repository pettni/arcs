function [V, K] = win_until(ts, B, P, quant1)
    % Compute the winning set of
    %  B U P
    % under (quant1, forall)-controllability
    V = uint32([]);
    K = containers.Map('KeyType', 'uint32', 'ValueType', 'any');
    while true
        [preV, preK] = ts.pre(V, 1:ts.n_a, quant1, 'forall');
        K = [preK; K];   % order important, K takes priority!
        Vt = union(P, intersect(B, preV));
        Vt = reshape(Vt, 1, length(Vt));
        for i=1:length(ts.pg_U)
            % Progress groups
            [preVinv, preKinv] = ts.pginv(ts.pg_U{i}, ts.pg_G{i}, V, B, quant1);
            K = [preKinv; K];   % order important, K takes priority!
            Vt = union(Vt, preVinv);
            Vt = reshape(Vt, 1, length(Vt));
        end
        if length(V) == length(Vt)
            break
        end
        V = Vt;
    end

    if ts.b_debug
        assert(all(ismember(setdiff(V, P), cell2mat(K.keys))))
    end
end
