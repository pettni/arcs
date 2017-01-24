function [ret, K] = pre(ts, X, U, quant1, quant2)
    % Compute pre(X) under (quant1, quant2)-controllability
    % and action set U
    K = containers.Map('KeyType', 'uint32', 'ValueType', 'any');

    pre_all = uint32([]);

    if ts.fast_enabled
        for i=1:length(X)
            pre_all = [pre_all ts.fast_pre(X(i))];
        end
    else
        for i=1:ts.num_trans()
            if ismember(ts.state2(i), X) && ismember(ts.action(i), U)
                pre_all(end+1) = ts.state1(i);
            end
        end
    end
    pre_all = unique(pre_all); 

    ret = uint32([]);
    for i = 1:length(pre_all)
        q = pre_all(i);     % candidate state
        act_list = zeros(1, ts.n_a);   % outcome per action 
        for a = 1:ts.n_a
            aPost = ts.post(q, a);
            if strcmp(quant2, 'exists')
                act_list(a) = any(ismember(aPost, X));
            else
                act_list(a) = all(ismember(aPost, X)) && length(aPost) > 0;
            end
        end
        if strcmp(quant1, 'exists')
            if any(act_list)
                ret(end+1) = q;
                K(q) = find(act_list); 
            end
        else
            if all(act_list)
                ret(end+1) = q;
                K(q) = 1:ts.n_a;
            end
        end
    end

    if ts.b_debug
        assert(all(ismember(ret, cell2mat(K.keys))))
    end
end
