function array = ts2array(ts)
% build the state-action transition matrix
    n_a = double(ts.n_a);
    n_s = double(ts.n_s);
    array = cell(ts.n_a,1);
    for i = 1:n_a
        idx = ts.action==i;
        % change sub to idx s.t. fast to assign 1's
        idx_m = sub2ind([n_s,n_s],ts.state1(idx),ts.state2(idx));
        array{i} = sparse([],[],[],n_s,n_s,length(idx_m));
        array{i}(idx_m)=true;
    end
end