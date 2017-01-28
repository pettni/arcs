function [W, K] = pre_pg(ts, V, B, quant1)
  % pre_pg: pre(V) under (quant1, forall) while remaining in B using progress groups
  % 
  % Returns a sorted set
  W = uint32([]);
  K = containers.Map('KeyType', 'uint32', 'ValueType', 'any');
  for i=1:length(ts.pg_U)
    % Progress groups
    if nargout > 1
      [preVinv, preKinv] = ts.pginv(ts.pg_U{i}, ts.pg_G{i}, V, B, quant1);
      K = [preKinv; K];   % order important, K takes priority!
    else
      preVinv = ts.pginv(ts.pg_U{i}, ts.pg_G{i}, V, B, quant1);
    end
    W = union(W, preVinv);
    W = reshape(W, 1, length(W));
  end
end