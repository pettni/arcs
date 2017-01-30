function [V, Cv, cont] = win_until(ts, B, P, quant1)
  % Compute the winning set of
  %  B U P
  % under (quant1, forall)-controllability
  %
  % Returns a sorted set
  %
  % Exanding algo

  V = uint32([]);
  Vlist = {};
  Klist = {};
  while true
    if nargout > 2
      C_rec = [];
      % Normal pre
      [preV, preK] = ts.pre(V, 'all', quant1, 'forall');
      Vlist{end+1} = preV;
      Klist{end+1} = preK;
      Vt = union(P, intersect(B, preV));

      % PG pre
      [preVinv, CpreVinv, preKinv] = ts.pre_pg(Vt, B, quant1);
      if ~isempty(setdiff(preVinv, Vt))
        Vlist{end+1} = preVinv;
        Klist{end+1} = preKinv;
      end
      C_rec = union(C_rec, CpreVinv);
      Vt = union(Vt, preVinv);
    elseif nargout > 1
      C_rec = [];
      % Normal pre
      preV = ts.pre(V, 'all', quant1, 'forall');
      Vt = union(P, intersect(B, preV));

      % PG pre
      [preVinv, CpreVinv] = ts.pre_pg(Vt, B, quant1);
      C_rec = union(C_rec, CpreVinv);
      Vt = union(Vt, preVinv);
    else
      preV = ts.pre(V, 'all', quant1, 'forall');
      Vt = union(P, intersect(B, preV));
      preVinv = ts.pre_pg(V, B, quant1);
      Vt = union(Vt, preVinv);
    end

    Vt = reshape(Vt, 1, length(Vt));

    if length(V) == length(Vt)
      break
    end

    V = Vt;
  end

  if nargout > 1
    Cv = union(C_rec, setdiff(ts.pre(V, 'all', quant1, 'forall'), V));
  end

  if nargout > 2
    cont = Controller(Vlist, Klist, 'reach', 'win_until');
  end
end
