function [W, cont] = pginv(ts, U, G, Z, B, quant1)
  % Compute U-controlled set
  % contained in G \cap B \setdiff Z that can be
  % used to force a transition to Z using the progress
  % group (U,G) under (quant1, forall)-controllability
  %  
  % Returns a sorted set

  W = setdiff(intersect(uint32(G), uint32(B)), uint32(Z));

  if ts.b_disable_pg || ...  % pgs disabled
     (strcmp(quant1, 'forall') && ~isempty(setdiff(1:ts.n_a, U))) || ... % uncontrolled actions
     isempty(intersect(ts.pre(Z, U, quant1, 'exists'), W))  % no reach to Z
    W = [];
    cont = Controller(W, containers.Map(), 'simple');
    return
  end

  while true
    if nargout > 1
      % union is sorted
      [preW, K] = ts.pre(union(W, Z), U, quant1, 'forall');
    else
      % union is sorted
      preW = ts.pre(union(W, Z), U, quant1, 'forall');
    end
    Wt = intersect(W, preW);
    if length(W) == length(Wt)
      break
    end
    W = Wt;
  end

  if nargout > 1
    [~, cont] = ts.pre(union(W, Z), U, quant1, 'forall');
    cont.restrict_to(W);
    cont.from = 'pginv';
  end
end
