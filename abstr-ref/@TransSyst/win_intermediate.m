function [V, Cv, cont] = win_intermediate(ts, A, B, P, C_list, quant1)
  % Compute winning set of
  %  []A && ( (B U P) || [] (B &&_i <>C_i) )
  % under (quant1, forall)-controllability
  %
  % Note: A must be sorted
  % Returns a sorted set
  % 
  % Contracting algorithm
  
  if strcmp(ts.sys_setting, TransSyst.bdd_set)
    if nargout == 1
      V = ts.bdd_sys.win_intermediate(A, B, P, C_list, quant1);
    elseif nargout == 2
      [V, Cv] = ts.bdd_sys.win_intermediate(A, B, P, C_list, quant1);
    elseif nargout == 3
      [V, Cv, cont] = ts.bdd_sys.win_intermediate(A, B, P, C_list, quant1);
    end
    return;
  end

  V = uint32(1:ts.n_s);
  Klist = {};

  iter = 1;
  while true
    Vt = V;
    preV = ts.pre(V, [], quant1, false);
    
    if nargout > 1
      Cv = [];
    end

    for i=1:length(C_list)
      Qi = union(P, intersect(intersect(B, C_list{i}), preV));
      Qi = reshape(Qi, 1, length(Qi));
      if nargout > 2
        [Vti, Cvi, Ki] = ts.win_until_and_always(A, B, Qi, quant1);
        Klist{i} = Ki;
      elseif nargout > 1
        [Vti, Cvi] = ts.win_until_and_always(A, B, Qi, quant1);
        Cv = union(Cv, Cvi);
      else
        Vti = ts.win_until_and_always(A, B, Qi, quant1);
      end
      Vt = intersect(Vt, Vti);
    end
    
    if nargout > 1 && iter == 1
      V1 = Vt;
      C_rec = Cv;
    end

    if length(V) == length(Vt)
      break
    end
    V = Vt;
    iter = iter + 1;
  end

  if nargout > 1
    % Contracting: C_rec U (V_1\V_last)
    Cv = union(C_rec, setdiff(V1, V));
  end

  if nargout > 2
    Vlist = {V};
    for i = 1:length(C_list)
      Vlist{end+1} = intersect(B, C_list{i});
    end
    cont = Controller(Vlist, Klist, 'recurrence', 'win_intermediate');
  end
end
