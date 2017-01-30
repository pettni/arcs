% Should be replaced when dual algos implemented!
function [V, cont] = win_always(ts, V, quant1)
  % Compute the winning set of
  %  [] V
  % under (quant1, forall)-controllability

  % Note: V must be sorted
  % Returns a sorted set
  while true
    Vt = intersect(V, ts.pre(V, 1:ts.n_a, quant1, 'forall'));
    if length(V) == length(Vt)
      break
    end
    V = Vt;
  end

  if nargout > 1
    [~, cont] = ts.pre(V, 1:ts.n_a, quant1, 'forall');
    cont.from = 'win_always';
  end

end
