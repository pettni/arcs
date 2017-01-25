function [trans trans_out] = compute_transitions_nld(fx, part, vars, dBound)
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
            [trans(:,:,i) trans_out(:,:,i)] = compute_transitions_nld(fx{i}, part, vars, dBound);
        end
        return;
    end

    N = length(part);
    trans = sparse(N, N);
    trans_out = zeros(1,N);
    for i=1:N
        [adj, dim] = part.get_neighbors(i);
        for j=adj
            if isTransNLinD(part(i), part(j), fx, vars, dBound)
                trans(i,j) = 1;
            end
        end
        if isTransOutNLinD(part(i), part.domain, fx, vars, dBound)
            trans_out(i) = 1;
        end
        trans(i,i) = not(isTransientNLinD(part(i),fx, vars, deg, dBound));
    end
end