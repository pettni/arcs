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
  p = Polynomial(p_coef, p_mon);

  g_mon =  [0  2  0;
            0  0  2];
  g_coef = [1 -1 -1];
  g = Polynomial(g_coef, g_mon);

  % Total equation is in 2 vars
  A_g = PolyLinTrans.eye(1, 2, 0, tot_deg).as_vector_trans;
  A_s = PolyLinTrans.mul_pol(2, sigma_deg, g).as_matrix_trans;
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

  % Constraint dimensions
  m_1 = size(A_g, 1);

  prob.a = [A_g             -A_s              -A_k];
  % These are equality constraints
  prob.blc  = [p.mon_vec(tot_deg)];  
  prob.buc  = [p.mon_vec(tot_deg)]; 

  % sdsos constraints
  prob = add_sdsos(prob, n_g+1, n_g+n_s);
  prob = add_sdsos(prob, n_g+n_s+1, n_g+n_s+n_k);

  % Objective
  prob.c = zeros(1, size(prob.a,2));
  prob.c(1) = 1;

  [r,res]=mosekopt('minimize echo(0)', prob); 

  gamma_opt = res.sol.itr.xx(1);

  verifyEqual(testCase, gamma_opt, 8.5, 'reltol', 0.05)

end