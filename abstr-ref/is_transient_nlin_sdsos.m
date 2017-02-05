function result = is_transient_nlin_sdsos(rec1, dyn_list, drec, tot_deg)
  % Determine if the modes dyn_list = {fx1, fx2, ...}, where
  % fxi is a vector of Polynomial, is transient on rec1 using 
  % relaxation order tot_deg
  %
  epsilon = 1;

  M = length(dyn_list);  % number of modes
  deg_f = 0;
  for m = 1:M
    deg_f = max([deg_f, dyn_list{m}.deg]);
  end

  if isempty(drec)
    % No disturbance
    drec = Rec(zeros(2,0));
  end

  all_rec = rec1 * drec;

  n_x = rec1.dim;
  n_d = drec.dim;

  if n_x+n_d ~= all_rec.dim
    error('dimension mismatch')
  end

  % For all m: -d B . f_m (x,d) - eps - \sum_i sji(x,d) gji(x,d)  >= 0 
  % Overall equation is in vars x,d and degree tot_deg

  tot_deg = tot_deg - mod(tot_deg, 2);  % make it even
  deg_sigma = tot_deg - 2;  % degree of all g's is 2
  deg_sigma = deg_sigma - mod(deg_sigma, 2);  % make it even
  deg_B = tot_deg + 1 - deg_f;

  % Variable vector  [B    K1     sigma1  K2     sigma2 ... Km  sigmam]
  %                  (vec) (smat) (smat)  (smat)
  
  n_b = count_monomials_leq(n_x, deg_B);  % # variables in B
  m_b = count_monomials_leq(n_x + n_d, tot_deg);  % # constraints in one block

  % Id transformation from K (smat)
  A_k = PolyLinTrans.eye(n_x+n_d, n_x+n_d, tot_deg, tot_deg).as_matrix_trans;
  n_k = size(A_k, 2);

  % Multiplication with g's transformation from sigma's (smat)
  A_si = cell(1, all_rec.dim);
  n_si = zeros(1, all_rec.dim);
  for i = 1:all_rec.dim
    % Create polynomial x_i - x_i^2
    mons = zeros(all_rec.dim,2);
    mons(i,:) = [1 2];
    p = Polynomial([1 -1], mons);
    A_si{i} = PolyLinTrans.mul_pol(...
                all_rec.dim, deg_sigma, ...
                tot_deg, p ).as_matrix_trans;
    n_si(i) = size(A_si{i}, 2);
  end
  A_s = [A_si{:}];
  n_s = sum(n_si);

  % Build Lie derivative operators B(x) -> dB.f_m (x,d)
  T_dB_list = cell(1,M);
  for m = 1:M
    T_dB = PolyLinTrans(n_x, n_x+n_d, deg_B, tot_deg);
    for j = 1:n_x
      fm_j = dyn_list{m}(j);

      % Get fm_j( x^- + (x^+ - x^-) x ) / (x_j^+ - x_j^-)
      fm_j_s = fm_j.shift(all_rec.xmin).scale(all_rec.xmax - all_rec.xmin) / ...
               (all_rec.xmax(j) - all_rec.xmin(j));

      T_mul = PolyLinTrans.mul_pol(n_x + n_d, deg_B-1, tot_deg, fm_j_s);
      T_lift = PolyLinTrans.eye(n_x, n_x + n_d, deg_B-1, deg_B-1);
      T_diff = PolyLinTrans.diff(n_x, deg_B, j);
      T_dB = T_dB + T_mul * T_lift * T_diff;
    end
    T_dB_list{m} = T_dB;
  end

  % Build matrix
  prob = [];
  prob.a = sparse(0,n_b);
  prob.blc = [];
  prob.buc = [];

  for m = 1:M
    A_dBm = T_dB_list{m}.as_vector_trans;

    % Stack  [   oldA     0        0;  = [oldb;
    %         -A_dBm  0   -A_k   -A_s]     eps]
    prob.a = [prob.a   zeros(size(prob.a, 1), n_k + n_s);
              -A_dBm zeros(m_b, size(prob.a, 2) - n_b) -A_k -A_s];
    prob.blc = [prob.blc; epsilon; zeros(m_b-1, 1)];
    prob.buc = [prob.buc; epsilon; zeros(m_b-1, 1)];
  end

  % Add sdsos constraints
  for i = 1:M
    idx1 = n_b + (n_k + n_s) * (i-1) + n_k;
    prob = add_sdsos(prob, n_b + (n_k + n_s) * (i-1) + 1, ...
                     idx1);   % make K_m sdd
    for j = 1:all_rec.dim
      prob = add_sdsos(prob, idx1 + sum(n_si(1:j-1)) + 1, ...
                       idx1 + sum(n_si(1:j)));  % make sigma_m sdd
    end
  end

  [~, res] = mosekopt('minimize echo(0)', prob); 

  result = strcmp(res.sol.itr.prosta, 'PRIMAL_AND_DUAL_FEASIBLE');

end