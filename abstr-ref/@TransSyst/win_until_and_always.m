function [V, Cv, cont] = win_until_and_always(ts, A, B, P, quant1)
  % Compute the winning set of
  %   []A && B U P
  % under (quant1, forall)-controllability
  %
  % Note: A must be sorted
  % Returns a sorted set
  if ~isempty(setdiff(1:ts.n_s, A))
    % Need to worry about []A
    if nargout > 2
      [Vinv, Cvinv, cont_inv] = ts.win_intermediate(uint32(1:ts.n_s), A, [], {uint32(1:ts.n_s)}, quant1);
      [V, Cvu, cont] = ts.win_until(intersect(B, Vinv), intersect(P, Vinv), quant1);
      cont.from = 'win_until_and_always';
      Cv = union(Cvinv, Cvu);
    elseif nargout > 1
      [Vinv, Cvinv] = ts.win_intermediate(uint32(1:ts.n_s), A, [], {uint32(1:ts.n_s)}, quant1);
      [V, Cvu] = ts.win_until(intersect(B, Vinv), intersect(P, Vinv), quant1);
      Cv = union(Cvinv, Cvu);
    else
      Vinv = ts.win_intermediate(uint32(1:ts.n_s), A, [], {uint32(1:ts.n_s)}, quant1);
      V = ts.win_until(intersect(B, Vinv), intersect(P, Vinv), quant1);
    end
  else
    % No need to worry about []A
    if nargout > 2
      [V, Cv, cont] = ts.win_until(B, P, quant1);
      cont.from = 'win_until';
    elseif nargout > 1
      [V, Cv] = ts.win_until(B, P, quant1);
    else
      V = ts.win_until(B, P, quant1);
    end
  end
end
