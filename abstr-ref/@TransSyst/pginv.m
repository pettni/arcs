function [W, Cw, cont] = pginv(ts, U, G, Z, B, quant1)
  % Compute U-controlled set
  % contained in G \cap B \setdiff Z that can be
  % used to force a transition to Z using the progress
  % group (U,G) under (quant1, forall)-controllability
  %  
  % Returns a sorted set
  % 
  % Contracting algorithm

  W1 = setdiff(intersect(uint32(G), uint32(B)), uint32(Z));
  W = W1;

  if ts.b_disable_pg || ...  % pgs disabled
     (~quant1 && ~isempty(setdiff(1:ts.n_a, U))) || ... % uncontrolled actions
     isempty(intersect(ts.pre(Z, U, quant1, 1), W))  % no reach to Z
    W = [];
    Cw = [];
    cont = Controller(W, containers.Map(), 'simple');
    return
  end

  while true
    % union is sorted
    preW = ts.pre(union(W, Z), U, quant1, false);
    Wt = intersect(W, preW);
    if length(W) == length(Wt)
      break
    end
    W = Wt;
  end

  if nargout > 1
    % Contracting: first set minus final set
    Cw = setdiff(W1, W);
  end

  if nargout > 2
    [~, cont] = ts.pre(union(W, Z), U, quant1, false);
    cont.restrict_to(W);
    cont.from = 'pginv';
  end
end
