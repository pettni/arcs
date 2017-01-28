function [V, K] = win_until(ts, B, P, quant1)
  % Compute the winning set of
  %  B U V
  % under (quant1, forall)-controllability
  %
  % Returns a sorted set

  V = uint32([]);
  K = containers.Map('KeyType', 'uint32', 'ValueType', 'any');
  while true
    if nargout > 1
      [preV, preK] = ts.pre(V, 1:ts.n_a, quant1, 'forall');
      [preVinv, preKinv] = ts.pre_pg(V, B, quant1);
      K = [preKinv; preK; K];   % order important, K takes priority!
    else
      preV = ts.pre(V, 1:ts.n_a, quant1, 'forall');
      preVinv = ts.pre_pg(V, B, quant1);
    end
    Vt = union(union(P, intersect(B, preV)), preVinv);
    Vt = reshape(Vt, 1, length(Vt));

    if length(V) == length(Vt)
      break
    end
    V = Vt;
  end
end
