function [V, Klist] = win_intermediate(ts, A, B, P, C_list, quant1)
    % Compute winning set of
    %  []A && ( (B U P) || [] (B &&_i <>C_i) )
    % under (quant1, forall)-controllability
    %
    % Note: A must be sorted
    % Returns a sorted set
    V = uint32(1:ts.n_s);
    Klist = {};
    while true
        Vt = V;
        preV = ts.pre(V, 1:ts.n_a, quant1, 'forall');
        for i=1:length(C_list)
            Qi = union(P, intersect(intersect(B, C_list{i}), preV));
            Qi = reshape(Qi, 1, length(Qi));
            if nargout > 1
                [Vti, Ki] = ts.win_until_and_always(A, B, Qi, quant1);
                Klist{i} = Ki;
            else
                Vti = ts.win_until_and_always(A, B, Qi, quant1);
            end
            Vt = intersect(Vt, Vti);
        end
        
        if length(V) == length(Vt)
            break
        end
        V = Vt;
    end
end
