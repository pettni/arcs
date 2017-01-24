% Should be replaced when dual algos implemented!
function [V, K] = win_always(ts, B, quant1)
    % Compute the winning set of
    %  [] B
    % under (quant1, forall)-controllability
    V = uint32(B);
    while true
        [preV, K] = ts.pre(V, 1:ts.n_a, quant1, 'forall');
        Vt = intersect(V, preV);
        if length(V) == length(Vt)
            break
        end
        V = Vt;
    end
end
