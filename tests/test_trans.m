function tests = exampleTest
  tests = functiontests(localfunctions);
end

function test_trans1(testCase)
	r1 = Rec([0 0; 1 1]);
	r2 = Rec([1 0; 2 1]);
	A = [0.1 0; 0 0.1];
	K = [1; 0];
	verifyEqual(testCase, is_trans(r1, r2, {A, K}), true);
	verifyEqual(testCase, is_trans(r2, r1, {A, K}), false);

	x = sdpvar(2,1);
	verifyEqual(testCase, is_trans(r1, r2, {A*x + K, x}), true);
	verifyEqual(testCase, is_trans(r2, r1, {A*x + K, x}), false);
end

function test_trans12(testCase)
  r1 = Rec([-1 -1; -0.1 0]);
  r2 = Rec([-1 0; -0.1 1]);

  A = [-1 -1; 1 -1];
  K = [1; 0];
  E = [0; 1];
  drec = Rec([0 0.2]);

  verifyEqual(testCase, is_trans(r1, r2, {A, K}), false);
  verifyEqual(testCase, is_trans(r2, r1, {A, K, E, drec}), true);

  x = sdpvar(2,1);
  d = sdpvar(1,1);
  verifyEqual(testCase, is_trans(r1, r2, {A*x + K, x}), false);
  verifyEqual(testCase, is_trans(r2, r1, {A*x + E*d + K, x, d, drec}), true);
end

function test_nonl(testCase)

  sdpvar x y;

  f = [-x-y^2+1; x*y-y];
  r1 = Rec([-1 0; 0 0.8]);
  r2 = Rec([0 0; 1 0.8]);
  verifyEqual(testCase, is_trans(r1, r2, {f, [x; y]}), true);
  verifyEqual(testCase, is_trans(r2, r1, {f, [x; y]}), false);

  r1 = Rec([-1 0; 0 1.2]);
  r2 = Rec([0 0; 1 1.2]);
  verifyEqual(testCase, is_trans(r2, r1, {f, [x; y]}), true);

end