function tests = exampleTest
  tests = functiontests(localfunctions);
end

function test_sdd(testCase)
  tot_deg = 8;
  sigma_deg = tot_deg - 2;  % 2 is deg of g

  % make it even
  sigma_deg = sigma_deg - mod(sigma_deg, 2);

  p_coef = [1 2 3 4 5];
  p_mon =  [0 1 0 1 2;
            1 0 2 1 0];

  g_mon =  [0  2  0;
            0  0  2];
  g_coef = [1 -1 -1];

  % Total equation is in 2 vars
  A_g = PolyLinTrans.eye(1, 2, 0, tot_deg).as_vector_trans;
  A_s = PolyLinTrans.mul_pol(2, sigma_deg, g_mon, g_coef).as_matrix_trans;
  A_k = PolyLinTrans.eye(2, 2, tot_deg, tot_deg).as_matrix_trans;

  [r, res] = mosekopt('symbcon echo(0)'); 
  prob          = []; 

  % Three constraints:
  %  1. gamma - sigma * g - K = p
  %  2. sigma = Ms
  %  3. K = Mk
  % 
  % Variable vector [gamma  sigma  K  Ms  Mk]

  % [Ag  -As  -Ak    0   0;    [p;
  %  0    -I    0   Ds   0;  =  0
  %  0     0   -I    0  DK]     0]

  % Variable dimensions
  n_g = 1;
  n_s = size(A_s, 2);
  n_k = size(A_k, 2);

  D_s = sdd_mat(n_s);
  D_k = sdd_mat(n_k);

  n_ms = size(D_s, 2);
  n_mk = size(D_k, 2);

  % Constraint dimensions
  m_1 = size(A_g, 1);
  m_2 = size(D_s, 1);
  m_3 = size(D_k, 1);

  prob.a = [A_g             -A_s              -A_k              sparse(m_1, n_ms) sparse(m_1, n_mk);
            sparse(m_2, n_g) -speye(m_2)       sparse(m_2, n_k) D_s               sparse(m_2, n_mk);
            sparse(m_3, n_g) sparse(m_3, n_s) -speye(m_3)       sparse(m_3, n_ms) D_k];
  % These are equality constraints
  prob.blc  = zeros(m_1 + m_2 + m_3, 1);  
  prob.buc  = zeros(m_1 + m_2 + m_3, 1); 

  for i=1:length(p_coef)
    prob.blc(mono_rank_grlex(2, p_mon(:,i)), 1) = p_coef(i);
    prob.buc(mono_rank_grlex(2, p_mon(:,i)), 1) = p_coef(i);
  end

  num_cone = (n_ms + n_mk)/3;

  % sdsos constraints
  prob.cones.type = res.symbcon.MSK_CT_RQUAD * ones(1, num_cone);
  prob.cones.sub = (1:3*num_cone) + n_g + n_s + n_k;
  prob.cones.subptr = 1:3:3*num_cone;

  % Objective
  prob.c = zeros(1, n_g+n_s+n_k+n_ms+n_mk);
  prob.c(1) = 1;

  % Unbounded vars
  prob.blx = -inf*ones(1, n_g+n_s+n_k+n_ms+n_mk);
  prob.bux = inf*ones(1, n_g+n_s+n_k+n_ms+n_mk);

  % param.MSK_IPAR_INTPNT_BASIS   = sc.MSK_OFF

  [r,res]=mosekopt('minimize echo(0)', prob); 

  gamma_opt = res.sol.itr.xx(1);

  verifyEqual(testCase, gamma_opt, 8.5, 'reltol', 0.05)

end