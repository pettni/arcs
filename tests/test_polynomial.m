function tests = exampleTest
  tests = functiontests(localfunctions);
end

function test_scale(testCase)
  sdpvar x y z;
  p = Polynomial(x^2 + x*y*z, [x; y; z]);

  q = p.scale([3 2 2]');

  verifyEqual(testCase, q, Polynomial(9*x^2 + 12*x*y*z, [x; y; z]))

end

function test_shift(testCase)
  sdpvar x y z;
  p = Polynomial(x^2 + x*y*z, [x; y; z]);

  q = p.shift([3 2 2]');
  q.reduce;

  verifyEqual(testCase, q, Polynomial((x+3)^2 + (x+3)*(y+2)*(z+2), [x; y; z]))

  p = Polynomial(1 + x^2 + x*y*z - z^3, [x; y; z]);
  
  q = p.shift([-2 2 2]');
  q.reduce;

  verifyEqual(testCase, q, Polynomial(1 + (x-2)^2 + (x-2)*(y+2)*(z+2) - (z+2)^3, [x; y; z]))

  p = Polynomial((x-5)^3 + (y-6)^8 + (z-12)^12, [x; y; z]);
  q = p.shift([5; 6; 12]);
  q.reduce;

  verifyEqual(testCase, q, Polynomial(x^3 + y^8 + z^12, [x;y;z]));
end