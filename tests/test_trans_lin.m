function tests = exampleTest
  tests = functiontests(localfunctions);
end

function test_trans1(testCase)
    r1 = Rec([0 0; 1 1]);
    r2 = Rec([1 0; 2 1]);
    A = [0.09 0; 0 0.09];
    K = [1; 0];
    verifyEqual(testCase, is_trans_lin(r1, r2, {A, K}), true);
    verifyEqual(testCase, is_trans_lin(r2, r1, {A, K}), false);
end

function test_trans12(testCase)
  r1 = Rec([-1 -1; -0.1 0]);
  r2 = Rec([-1 0; -0.1 1]);

  A = [-1 -1; 1 -1];
  K = [1; 0];
  E = [0; 1];
  drec = Rec([0 0.2]);

  verifyEqual(testCase, is_trans_lin(r1, r2, {A, K}), false);
  verifyEqual(testCase, is_trans_lin(r2, r1, {A, K, E}, drec), true);
end

function setupOnce(testCase)  % do not change function name
  global ops
  ops = sdpsettings('solver', 'mosek', 'cachesolvers', 1, 'verbose', 0);
end
