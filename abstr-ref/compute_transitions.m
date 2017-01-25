function [trans trans_out] = compute_transitions(fx, part)
    % Compute transition matrix for the partition part
    % Returns a sparse matrix M where M(i,j)=1 if there
    % is a transition from cell i to cell j.
    %

    if length(fx)>1
        % If there are several modes, compute transitions for each one
        trans = zeros(length(part), length(part), length(fx));
        trans_out = zeros(length(part), 1, length(fx));
        for i=1:length(fx)
            [trans(:,:,i) trans_out(:,:,i)] = compute_transitions(fx{i}, part);
        end
        return;
    end

    N = length(part);
    trans = sparse(N, N);
    trans_out = zeros(1,N);
    for i=1:N
        [adj, dim] = part.get_neighbors(i);
        for j=adj
            if isTransLin(part(i), part(j), fx)
                trans(i,j) = 1;
            end
        end
        if isTransOutLin(part(i), part.domain, fx)
            trans_out(i) = 1;
        end
        trans(i,i) = not(isTransientLin(part(i),fx));
    end
end