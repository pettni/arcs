%% Main function to generate tests
function tests = exampleTest
  tests = functiontests(localfunctions);
end

%% Test Functions
function test_transient1(testCase)
  r = Rec([0 0; 1 1]);
  A1 = [0 0; 0 0];
  K1 = [1; 0];
  A2 = [0 0; 0 0];
  K2 = [0; 1];

  dyn_list = {{A1, K1}, {A2, K2}};

  verifyEqual(testCase, is_transient_lin(r, dyn_list{1}), true);
  verifyEqual(testCase, is_transient_lin(r, dyn_list{2}), true);
end

function test_transient2(testCase)
  r = Rec([0 0; 1 1]);
  A1 = [0 0; 0 0];
  K1 = [1; 0];
  A2 = [0 0; 0 0];
  K2 = [-1; 0];

  dyn_list = {{A1, K1}, {A2, K2}};

  verifyEqual(testCase, is_transient_lin(r, dyn_list{1}), true);
  verifyEqual(testCase, is_transient_lin(r, dyn_list{2}), true);
end

function test_transient3(testCase)
  r = Rec([-1 -1 -1; 1 1 1]);
  A1 = [-1 0 0; 0 -1 0; 0 0 -1];
  K1 = [0; 0; 0.1];
  A2 = [-1 0 0; 0 -1 0; 0 0 -1];
  K2 = [0; 0; 0.1];

  dyn_list = {{A1, K1}, {A2, K2}};

  verifyEqual(testCase, is_transient_lin(r, dyn_list{1}), false);
  verifyEqual(testCase, is_transient_lin(r, dyn_list{2}), false);
end

function test_radiant1(testCase)
  r = Rec([20 28; 20 28; 20 28]);

  % These are from the CDC'14 paper
  % Crucial to divide by 10!
  A1 = [-0.0089    0.0020    0.0019;
        0.0040   -0.0073    0.0030;
        0.0040    0.0020   -0.0062];
  K1 = [0.0900;
    0.0107;
    0.1020/10];

  A2 = [-0.0039    0.0020    0.0019;
    0.0040   -0.0073    0.0030;
    0.0040    0.0020   -0.0062];
  K2 = [0;
    0.0107;
    0.1020/10];

  verifyEqual(testCase, is_transient_lin(r, {A1, K1}), true);
  verifyEqual(testCase, is_transient_lin(r, {A2, K2}), true);
end

function test_radiant2(testCase)
  r = Rec([20 28; 20 28; 20 28]);

  % These are the "newer ones"
  A1 = 1e-3 * [   -0.0413    0.0106    0.0080;
    0.4377   -0.4869    0.0260;
    0.4377    0.0346   -0.4955];
  K1 = [0.0004;
    0.0010;
    0.0011];

  A2 = 1e-3 * [-0.0186    0.0106    0.0080;
    0.4377   -0.4869    0.0260;
    0.4377    0.0346   -0.4955];
  K2 = [0;
    0.0010;
    0.0011];

  verifyEqual(testCase, is_transient_lin(r, {A1, K1}), true);
  verifyEqual(testCase, is_transient_lin(r, {A2, K2}), true);
end


function test_disturbance(testCase)
  r = Rec([-1 -1; 1 1]);
  A = [0 0; 0 0];
  B = [1; 0];
  E = [1; 0];

  drec = Rec([-1.1, 1.1]);

  verifyEqual(testCase, is_transient_lin(r, {A, B}), true);
  verifyEqual(testCase, is_transient_lin(r, {A, B, E, drec}), false);

end

function setupOnce(testCase)  % do not change function name
  global ops
  ops = sdpsettings('solver', 'mosek', 'cachesolvers', 1, 'verbose', 0);
end
