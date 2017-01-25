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
      K = [preK; K];   % order important, K takes priority!
    else
      preV = ts.pre(V, 1:ts.n_a, quant1, 'forall');
    end
    Vt = union(P, intersect(B, preV));
    Vt = reshape(Vt, 1, length(Vt));
    for i=1:length(ts.pg_U)
      % Progress groups
      if nargout > 1
        [preVinv, preKinv] = ts.pginv(ts.pg_U{i}, ts.pg_G{i}, V, B, quant1);
        K = [preKinv; K];   % order important, K takes priority!
      else
        preVinv = ts.pginv(ts.pg_U{i}, ts.pg_G{i}, V, B, quant1);
      end
      Vt = union(Vt, preVinv);
      Vt = reshape(Vt, 1, length(Vt));
    end
    if length(V) == length(Vt)
      break
    end
    V = Vt;
  end
end
