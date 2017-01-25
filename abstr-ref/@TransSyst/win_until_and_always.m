function [ret, K] = win_until_and_always(ts, A, B, P, quant1)
    % Compute the winning set of
    %   []A && B U P
    % under (quant1, forall)-controllability
    %
    % Note: A must be sorted
    % Returns a sorted set
    if nargout > 1
        [Vinv, Kinv] = ts.win_always(A, quant1);
        [ret, K] = ts.win_until(intersect(B, Vinv), intersect(P, Vinv), quant1);
        for i=1:length(ret)
            if ~ismember(ret(i), P)
                K(ret(i)) = intersect(K(ret(i)), Kinv(ret(i)));
            end
        end
    else
        Vinv = ts.win_always(A, quant1);
        ret = ts.win_until(intersect(B, Vinv), intersect(P, Vinv), quant1);
    end
end
