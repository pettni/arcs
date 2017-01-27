%% Main function to generate tests
function tests = exampleTest
  tests = functiontests(localfunctions);
end

function test_mldiv(testCase)
  r1 = Rec([0 0; 1 1]);
  r2 = Rec([0.5 0.5; 2 2]);

  d = mldivide(r1,r2);

  verifyEqual(testCase, length(d), 2)
  verifyEqual(testCase, isInside(d, [0.1; 0.1]), true)
  verifyEqual(testCase, isInside(d, [0.9; 0.9]), false)

  r3 = Rec([-1 -1; 0.1 0.1]);
  d = mldivide(d, r3);

  verifyEqual(testCase, length(d), 3)
  verifyEqual(testCase, isInside(d, [0.15; 0.15]), true)
  verifyEqual(testCase, isInside(d, [0.05; 0.05]), false)

  r4 = Rec([0 0; 1 1]);
  r5 = Rec([-1 0.5; 2 2]);
  r6 = Rec([0 0; 1 0.5]);
  verifyEqual(testCase, mldivide(r4, r5) == r6, true)
end

function test_isect(testCase)
  r1 = Rec([0 0; 1 1]);
  r2 = Rec([0.5 0.5; 2 2]);
  r3 = Rec([1 1; 2 2]);
  r4 = Rec([1.1 1.1; 2 2]);

  r5 = Rec([3 3; 4 4]);

  verifyEqual(testCase, intersects(r1, r2), true)
  verifyEqual(testCase, intersects(r1, r3), true)
  verifyEqual(testCase, intersects(r1, r4), false)

  verifyEqual(testCase, intersects([r1,r2], r4), true)
  verifyEqual(testCase, intersects(r4, [r1,r2]), true)
  verifyEqual(testCase, intersects([r1,r2], r5), false)
  verifyEqual(testCase, intersects(r5, [r1,r2]), false)
  verifyEqual(testCase, intersects([r1, r2], [r3, r4]), true)
  verifyEqual(testCase, intersects([r3, r4], [r1, r2]), true)
end