function [V, CV, Vlist, Klist] = win_primal(ts, A, B, C_list, quant1, V)
  % Compute winning set of
  %  []A && <>[]B &&_i []<>C_i
  % under (quant1, forall)-controllability
  % with the initial condition V ("warm start")
  %
  % Returns a sorted set
  
  if nargin<6
    V = [];
  end

  if isempty(A)
    A = uint32(1:ts.n_s);
  end
  if isempty(B)
    B = uint32(1:ts.n_s);
  end
  if isempty(C_list)
    C_list = {uint32(1:ts.n_s)};
  end

  Vlist = {};
  Klist = {};

  V = uint32(V);
  A = sort(A);

  ts.create_fast();

  while true
    Z = ts.pre(V, 1:ts.n_a, quant1, 'forall');
    for i=1:length(ts.pg_U)
      Z = union(Z, ts.pginv(ts.pg_U{i}, ts.pg_G{i}, V, A, quant1));
    end

    if nargout > 2
      [Vt, Kt] = ts.win_intermediate(A, B, Z, C_list, quant1);
    else
      Vt = ts.win_intermediate(A, B, Z, C_list, quant1);
    end

    if length(Vt) == length(V)
      break
    end

    if nargout > 2
      Klist{end+1} = Kt;
      Vlist{end+1} = Vt;
    end

    V = Vt;
  end

  % Candidate set
  CV = setdiff(ts.pre(V, 1:ts.n_a, quant1, 'exists'), V);

  CV = union(CV, setdiff(intersect(A, B), ts.win_always(intersect(A, B), quant1)));

  % for i=1:length(C_list)
  %   W1j = ts.win_until(B, intersect(B, C_list{i}), quant1);
  %   CV = union(CV, W1j);
  % end

  % CV = union(CV, setdiff(A, ts.win_always(A, quant1)));

end
