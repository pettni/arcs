function [V, Cv, cont] = win_until_and_always(ts, A, B, P, quant1)
  % Compute the winning set of
  %   []A && B U P
  % under (quant1, forall)-controllability
  %
  % Note: A must be sorted
  % Returns a sorted set
  if nargout > 2
    [Vinv, Cvinv, cont_inv] = ts.win_always(A, quant1);
    [V, Cvu, cont] = ts.win_until(intersect(B, Vinv), intersect(P, Vinv), quant1);
    cont.add_safety(cont_inv);
    cont.from = 'win_until_and_always';
  end

  if nargout > 1
    [Vinv, Cvinv] = ts.win_always(A, quant1);
    [V, Cvu] = ts.win_until(intersect(B, Vinv), intersect(P, Vinv), quant1);
    Cv = union(Cvinv, Cvu);
  else
    Vinv = ts.win_always(A, quant1);
    V = ts.win_until(intersect(B, Vinv), intersect(P, Vinv), quant1);
  end
end
