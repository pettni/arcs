function [W, Cw, cont] = pre_pg(ts, V, B, quant1)
  % pre_pg: pre(V) under (quant1, forall) while remaining in B using progress groups
  % 
  % Returns a sorted set
  
  if strcmp(ts.sys_setting, TransSyst.bdd_set)
    if nargout == 1
      W = ts.bdd_sys.pre_pg(V, B, quant1);
    elseif nargout == 2
      [W, Cw] = ts.bdd_sys.pre_pg(V, B, quant1);
    elseif nargout == 3
      [W, Cw, cont] = ts.bdd_sys.pre_pg(V, B, quant1);
    end
    return;
  end
  
  W = uint32(V);
  Vlist = {V};
  Cw = [];
  Klist = {Controller(W, containers.Map(), 'simple')};

  for i=1:length(ts.pg_U)
    % Progress groups
    if nargout > 2
      [preVinv, Cw, preKinv] = ts.pginv(ts.pg_U{i}, ts.pg_G{i}, W, B, quant1);
      if length(preVinv) > 0
        Vlist{end+1} = preVinv;
        Klist{end+1} = preKinv;
      end
    elseif nargout > 1
      [preVinv, Cw] = ts.pginv(ts.pg_U{i}, ts.pg_G{i}, W, B, quant1);
    else
      preVinv = ts.pginv(ts.pg_U{i}, ts.pg_G{i}, W, B, quant1);
    end
    W = union(W, preVinv);
    W = reshape(W, 1, length(W));
  end

  if nargout > 2
    cont = Controller(Vlist, Klist, 'reach', 'pre_pg');
  end
end