function [V, Cv, cont] = win_dual(ts, A, B_list, C, quant1, quant2)
  % win_until_or_always: compute winning set of      
  %   <>A || or_i <>[]B_i || []<>C
  % under (quant1, quant2)-controllability

  if isa(quant1, 'char') && strcmp(quant1, 'exists')
    quant1_bool = true;
  elseif isa(quant1, 'char') && strcmp(quant1, 'forall')
    quant1_bool = false;
  else
    error('quantifier must be exists or forall')
  end

  if isa(quant2, 'char') && strcmp(quant2, 'exists')
    quant2_bool = true;
  elseif isa(quant2, 'char') && strcmp(quant2, 'forall')
    quant2_bool = false;
  else
    error('quantifier must be exists or forall')
  end

  if quant2_bool
    % quant2 is exists -- dualize
    cA = setdiff(uint32(1:ts.n_s), A);
    cC = setdiff(uint32(1:ts.n_s), C);
    cB_list = {};
    for i=1:length(B_list)
      cB_list{i} = setdiff(uint32(1:ts.n_s), B_list{i});
    end
    V = setdiff(uint32(1:ts.n_s), win_dual(cA, cB_list, cC, ~quant1_bool, ~quant2_bool));
    return
  end

  % quant2 is forall
  V = uint32(1:ts.n_s);

  while true
    preV = ts.pre(V, [], quant1_bool, quant2_bool);
    P_inner = union(A, intersect(C, preV));
    Vt = ts.win_eventually_or_persistence(P_inner, B_list, quant1_bool);
    if length(Vt) == length(V)
      break
    end
    V = Vt;
  end

  if nargout > 2
    [~, ~, cont] = ts.win_eventually_or_persistence(P_inner, B_list, quant1_bool);
    cont.from = 'win_dual';
  end
end