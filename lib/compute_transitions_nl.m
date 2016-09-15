function [trans trans_out] = compute_transitions_nl(fx, part, vars)
    % Compute transition matrix for the partition part
    % Returns a sparse matrix M where M(i,j)=1 if there
    % is a transition from cell i to cell j.
    %
    deg = 2;
    if length(fx)>1
        % If there are several modes, compute transitions for each one
        trans = zeros(length(part), length(part), length(fx));
        trans_out = zeros(length(part), 1, length(fx));
        for i=1:length(fx)
            [trans(:,:,i) trans_out(:,:,i)] = compute_transitions_nl(fx{i}, part, vars);
        end
        return;
    end

    N = length(part);
    trans = sparse(N, N);
    trans_out = zeros(1,N);
    for i=1:N
        [adj, dim] = part.get_neighbors(i);
        for j=adj
            if isTransNLin(part(i), part(j), fx, vars)
                trans(i,j) = 1;
            end
        end
        if isTransOutNLin(part(i), part.domain, fx, vars)
            trans_out(i) = 1;
        end
        trans(i,i) = not(isTransientNLin(part(i),fx, vars, deg));
    end
end