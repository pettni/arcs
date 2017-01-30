function [V, cont] = win_until(ts, B, P, quant1)
  % Compute the winning set of
  %  B U V
  % under (quant1, forall)-controllability
  %
  % Returns a sorted set

  V = uint32([]);
  K = containers.Map('KeyType', 'uint32', 'ValueType', 'any');
  Vlist = {};
  Klist = {};
  while true
    if nargout > 1
      [preV, preK] = ts.pre(V, 1:ts.n_a, quant1, 'forall');
      Vlist{end+1} = preV;
      Klist{end+1} = preK;
      Vt = union(V, intersect(B, preV));
      [preVinv, preKinv] = ts.pre_pg(Vt, B, quant1);
      Vlist{end+1} = preVinv;
      Klist{end+1} = preKinv;
      Vt = union(Vt, preVinv);
    else
      preV = ts.pre(V, 1:ts.n_a, quant1, 'forall');
      preVinv = ts.pre_pg(V, B, quant1);
      Vt = union(union(P, intersect(B, preV)), preVinv);
    end

    Vt = reshape(Vt, 1, length(Vt));

    if length(V) == length(Vt)
      break
    end
    V = Vt;
  end

  if nargout > 1
    cont = Controller(Vlist, Klist, 'reach');
  end
end
