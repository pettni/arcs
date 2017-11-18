function [V, C, cont] = win_eventually_or_persistence(ts, P, B_list, quant1)
  % win_until_or_always: compute winning set of      
  %    <>P || or_i <>[] B_i
  % under (quant1, 'forall')-controllability

  V = [];

  if nargout > 2
    Klist = {};
    Vlist = {};
  end

  while true
    if nargout > 2
      [preV, K] = ts.pre(V, [], quant1, 0);
      P_inner = union(P, ts.pre(V, [], quant1, 0));
      Vlist{end} = P_inner;
      Klist{end} = K;
    else
      preV = ts.pre(V, [], quant1, 0);
      P_inner = union(P, ts.pre(V, [], quant1, 0));
    end
    
    if nargout > 2
      [preVinv, ~, Kinv] = ts.pre_pg(V, uint32(1:ts.n_s), quant1);
      P_inner = union(P_inner, preVinv);
      Vlist{end} = P_inner;
      Klist{end} = Kinv;
    else
      preVinv = ts.pre_pg(V, uint32(1:ts.n_s), quant1);
      P_inner = union(P_inner, preVinv);
    end
    
    if nargout > 2
      [Vt, ~, Kt] = ts.win_until_or_always(B_list, P_inner, quant1);
      Vlist{end} = Vt;
      Klist{end} = Kt;
    else
      Vt = ts.win_until_or_always(B_list, P_inner, quant1);
    end
    
    if length(Vt) == length(V)
      break
    end

    V = Vt;
  end

  if nargout > 2
    cont = Controller(Vlist, Klist, 'reach', 'win_eventually_or_persistence');
  end
end
