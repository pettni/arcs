% Should be replaced when dual algos implemented!
function [V, Cv, cont] = win_always(ts, B, quant1)
  % Compute the winning set of
  %  [] B
  % under (quant1, forall)-controllability

  % Note: V must be sorted
  % Returns a sorted set
  % 
  % Contracting algorithm

  V = B;
  while true
    Vt = intersect(B, ts.pre(V, 'all', quant1, 'forall'));
    if length(V) == length(Vt)
      break
    end
    V = Vt;
  end

  if nargout > 1
    % Contracting: first minus last
    Cv = setdiff(B, V);
  end

  if nargout > 2
    [~, cont] = ts.pre(V, 'all', quant1, 'forall');
    cont.from = 'win_always';
  end

end
