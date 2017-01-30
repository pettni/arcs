function [V, Cv, cont] = win_primal(ts, A, B, C_list, quant1, V)
  % Compute winning set of
  %  []A && <>[]B &&_i []<>C_i
  % under (quant1, forall)-controllability
  % with the initial condition V ("warm start")
  %
  % Returns a sorted set
  %
  % Expanding algo
  
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

  iter = 1;
  while true
    Z = ts.pre(V, 'all', quant1, 'forall');
    Z = union(Z, ts.pre_pg(V, A, quant1));

    if nargout > 2
      [Vt, Ct, Kt] = ts.win_intermediate(A, B, Z, C_list, quant1);
    elseif nargout > 1
      [Vt, Ct] = ts.win_intermediate(A, B, Z, C_list, quant1);
    else
      Vt = ts.win_intermediate(A, B, Z, C_list, quant1);
    end

    if nargout > 1 && iter == 1
      C_rec = Ct;
    end

    if length(Vt) == length(V)
      break
    end

    if nargout > 2
      Klist{end+1} = Kt;
      Vlist{end+1} = Vt;
    end

    V = Vt;
    iter = iter+1;
  end

  % Candidate set
  if nargout > 1
    % Expanding: C_rec + pre(V)\V
    Cv = union(setdiff(ts.pre(V, 'all', quant1, 'exists'), V), C_rec);
  end

  % Controller
  if nargout > 2
    cont = Controller(Vlist, Klist, 'reach', 'win_primal');
  end

end
