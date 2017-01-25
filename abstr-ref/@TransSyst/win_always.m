% Should be replaced when dual algos implemented!
function [V, K] = win_always(ts, V, quant1)
  % Compute the winning set of
  %  [] V
  % under (quant1, forall)-controllability

  % Note: V must be sorted
  % Returns a sorted set
  while true
    if nargout > 1
      [preV, K] = ts.pre(V, 1:ts.n_a, quant1, 'forall');
    else
      preV = ts.pre(V, 1:ts.n_a, quant1, 'forall');
    end
    Vt = intersect(V, preV);
    if length(V) == length(Vt)
      break
    end
    V = Vt;
  end
end
