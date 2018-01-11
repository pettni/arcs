function result = is_trans_nlin_sdsos(rec1, rec2, fx, drec, tot_deg)
  % is_trans_nlin(rec1, rec2, dyn): Return false if a certificate that 
  % guarantees the non-existence of a flow from rec1 to rec2 under 
  % nonlinear dynamics dyn is found, true otherwise
  % 
  % Inputs 
  %   - Rec rec1, rec2: sets
  %   - Polynomial fx: dynamics \dot x = fx(x,d), d \in drec
  %   - double tot_deg: relaxation order for moment hierarchy

  result = true;

  irec = intersect(rec1,rec2);
  % overlapping
  if irec.isFullDim
    warning('is_trans_nlin_sdsos: was called with overlapping Recs - returning true');
    result = true;
    return;
  end

  % are they adjacent?
  if length(irec.getFlatDims) ~= 1
    result = false;
    return;
  end

  if isempty(drec)
    drec = Rec(zeros(0,2));
  end

  % find "flat" dimension
  flatdim = irec.getFlatDims;
  if rec1.xmax(flatdim)==rec2.xmin(flatdim)
    flatdim_val = rec1.xmax(flatdim);
    obj = fx(flatdim);
  else
    flatdim_val = rec1.xmin(flatdim);
    obj = -fx(flatdim);
  end 

  all_rec = irec.elvar(flatdim) * drec;
  red_dim = all_rec.dim;

  % Reduce dimension by eliminating flatdim
  T_el = PolyLinTrans.elvar(red_dim + 1, obj.deg, flatdim, flatdim_val);
  obj_red = T_el * obj;

  % Set up optimization problem
  prob = []; 

  % Basic constraint:
  %  1. gamma - \sum_i si * gxi - K = obj
  % 
  % Variable vector [gamma  K  s1 ... si]

  sigma_deg = tot_deg - 2;   % degree of all g's is 2

  Ag = PolyLinTrans.eye(red_dim, red_dim, 0, tot_deg).as_vector_trans;
  Ak = PolyLinTrans.eye(red_dim, red_dim, tot_deg, tot_deg).as_matrix_trans;

  n_g = 1;
  n_k = size(Ak, 2);

  Asi = cell(1, all_rec.dim);
  nsi = zeros(1, all_rec.dim);
  for i = 1:all_rec.dim
    g = all_rec.bounding_polynomial(i);
    Asi{i} = PolyLinTrans.mul_pol(all_rec.dim, sigma_deg, tot_deg, g).as_matrix_trans;
    nsi(i) = size(Asi{i}, 2);
  end

  % Add [Ag   -Ak   -As1 -As2 ... -Asn] =   [obj];
  As = [Asi{:}];
  prob.a = [Ag  -Ak  -As];
  prob.blc = obj_red.mon_vec(tot_deg);
  prob.buc = obj_red.mon_vec(tot_deg);

  % Make K sdsos
  prob = add_sdsos(prob, 1 + n_g, n_g + n_k);

  % Make the si's sdsos
  for i=1:red_dim
    prob = add_sdsos(prob, 1 + n_g + n_k + sum(nsi(1:i-1)), ...
                           n_g + n_k + sum(nsi(1:i)));
  end

  % Set objective
  prob.c = zeros(1, size(prob.a, 2));
  prob.c(1) = 1;

  % Optimize
  [~, res] = mosekopt('minimize echo(0)', prob); 

  gamma_opt = res.sol.itr.xx(1);

  if strcmp(res.sol.itr.solsta, 'OPTIMAL') && gamma_opt < 0
    result = false;
  end
end