function tests = exampleTest
  tests = functiontests(localfunctions);
end

function test_trans1(testCase)
	r1 = Rec([0 0; 1 1]);
	r2 = Rec([1 0; 2 1]);
	A = [0.09 0; 0 0.09];
	K = [1; 0];

	x = sdpvar(2,1);
	verifyEqual(testCase, is_trans_nlin(r1, r2, {A*x + K, x, []}, 4), true);
	verifyEqual(testCase, is_trans_nlin(r2, r1, {A*x + K, x, []}, 4), false);
end

function test_trans12(testCase)
  r1 = Rec([-1 -1; -0.1 0]);
  r2 = Rec([-1 0; -0.1 1]);

  A = [-1 -1; 1 -1];
  K = [1; 0];
  E = [0; 1];
  drec = Rec([0 0.2]);

  x = sdpvar(2,1);
  d = sdpvar(1,1);
  verifyEqual(testCase, is_trans_nlin(r1, r2, {A*x + K, x, []}, 4), false);
  verifyEqual(testCase, is_trans_nlin(r2, r1, {A*x + E*d + K, [x; d], drec}, 4), true);
end

function test_nonl(testCase)

  sdpvar x y;

  f = [-x-y^2+1; x*y-y];
  r1 = Rec([-1 0; 0 0.8]);
  r2 = Rec([0 0; 1 0.8]);
  verifyEqual(testCase, is_trans_nlin(r1, r2, {f, [x; y], []}, 4), true);
  verifyEqual(testCase, is_trans_nlin(r2, r1, {f, [x; y], []}, 4), false);

  r1 = Rec([-1 0; 0 1.2]);
  r2 = Rec([0 0; 1 1.2]);
  verifyEqual(testCase, is_trans_nlin(r2, r1, {f, [x; y], []}, 4), true);

end

function setupOnce(testCase)  % do not change function name
  global ops
  ops = sdpsettings('solver', 'mosek', 'cachesolvers', 1, 'verbose', 0);
end
