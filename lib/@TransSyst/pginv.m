function [W, K] = pginv(ts, U, G, Z, B, quant1)
    % Compute U-controlled set
    % contained in G \cap B \setdiff Z that can be
    % used to force a transition to Z using the progress
    % group (U,G) under (quant1, forall)-controllability

    if strcmp(quant1, 'forall') && ~isempty(setdiff(1:ts.n_a, U))
        W = [];
        K = containers.Map('KeyType', 'uint32', 'ValueType', 'any');
        return
    end

    W = setdiff(intersect(uint32(G), uint32(B)), uint32(Z));

    while true
        [preW, K] = ts.pre(union(W, Z), U, quant1, 'forall');
        Wt = intersect(W, preW);
        if length(W) == length(Wt)
            break
        end
        W = Wt;
    end

    if ts.b_debug
        assert(all(ismember(W, cell2mat(K.keys))))
    end
end
