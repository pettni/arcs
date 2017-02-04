function result = is_transient_nlin(rec1, dyn_list, tot_deg)
  % Determine if the modes dyn_list = {fx1, fx2, ...}, where
  % fx = {f, xvar, (dvar, drec)} is transient on rec1 using 
  % relaxation order tot_deg
  global ops;

  epsilon = 10;

  all_rec = rec1;
  all_var = dyn_list{1}{2};
  n_x = length(all_var);
  deg_f = 0;

  M = length(dyn_list);  % number of modes

  for i = 1:M
    if length(dyn_list{i}) > 2
      % Disturbance
      all_var = [all_var; dyn_list{i}{3}];
      all_rec = all_rec * dyn_list{i}{4};
    end
    deg_f = max(deg_f, degree(dyn_list{i}{1}));
  end

  n_d = length(all_var) - n_x;

  % For all j: -d B . f_j (x,d) - eps - \sum_i sji gji  >= 0 
  %                                            forall x,d \in all_rec 
  % Overall equation is in vars x,d and degree tot_deg

  tot_deg = tot_deg - mod(tot_deg, 2);  % should be even
  deg_sigma = tot_deg - 2;
  deg_sigma = deg_sigma - mod(deg_sigma, 2);  % should be even
  deg_B = tot_deg + 1 - deg_f;

  % Variable vector  [B    K1     sigma1  K2     sigma2 ... Km  sigmam]
  %                  (vec) (mat)  (mat)   (mat)
  n_b = count_monomials_leq(n_x, deg_B);

  m_b = count_monomials_leq(n_x + n_d, tot_deg);

  prob = [];
  prob.a = sparse(0,n_b);
  prob.blc = [];
  prob.buc = [];

  % Id transformation from K (smat)
  A_k = PolyLinTrans.eye(n_x+n_d, n_x+n_d, tot_deg, tot_deg).as_matrix_trans;
  n_k = size(A_k, 2);

  % Mul transformation from sigma (smat)
  A_si = cell(1, all_rec.dim);
  n_si = zeros(1, all_rec.dim);
  for i = 1:all_rec.dim
    g = all_rec.bounding_polynomial(i);
    A_si{i} = PolyLinTrans.mul_pol(all_rec.dim, deg_sigma, tot_deg, g).as_matrix_trans;
    n_si(i) = size(A_si{i}, 2);
  end
  A_s = [A_si{:}];
  n_s = sum(n_si);

  for m = 1:M
    % Build Lie derivative operator
    T_dB = PolyLinTrans(n_x, n_x+n_d, deg_B, tot_deg);
    for j = 1:n_x
      fm_j = Polynomial(dyn_list{m}{1}(j), all_var);
      T_mul = PolyLinTrans.mul_pol(n_x + n_d, deg_B-1, tot_deg, fm_j);
      T_lift = PolyLinTrans.eye(n_x, n_x + n_d, deg_B-1, deg_B-1);
      T_diff = PolyLinTrans.diff(n_x, deg_B, j);
      T_dB = T_dB + T_mul * T_lift * T_diff;
    end

    A_dB = T_dB.as_vector_trans;

    assert(size(A_dB, 2) == n_b);
    assert(size(A_dB, 1) == m_b);

    numvar = size(prob.a, 2);
    numcon = size(prob.a, 1);

    % Stack  [   oldA     0        0   = [oldb;
    %         -A_db  0   -A_k   -A_s]     eps]

    prob.a = [prob.a   zeros(size(prob.a, 1), n_k + n_s);
              -A_dB zeros(m_b, numvar - n_b) -A_k -A_s];
    prob.blc = [prob.blc; epsilon; zeros(m_b-1, 1)];
    prob.buc = [prob.buc; epsilon; zeros(m_b-1, 1)];
  end

  % Add sdsos constraints
  for i = 1:M
    idx0 = n_b + (n_k + n_s) * (i-1) + 1;
    idx1 = n_b + (n_k + n_s) * (i-1) + n_k;
    prob = add_sdsos(prob, idx0, idx1);   % make K_m sdd
    for j = 1:all_rec.dim
      idx2 = idx1 + sum(n_si(1:j-1)) + 1;
      idx3 = idx1 + sum(n_si(1:j));
      prob = add_sdsos(prob, idx2, idx3);  % make sigma_m sdd
    end
  end

  [r, res] = mosekopt('minimize echo(0)', prob); 

  % B = Polynomial(n_x, res.sol.itr.xx(1:n_b));
  % dB = T_dB * B;   % should be negative for all modes

  if strcmp(res.sol.itr.prosta, 'PRIMAL_AND_DUAL_FEASIBLE')
    result = true;
  else
    result = false;
  end

end