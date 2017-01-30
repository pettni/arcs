function tests = exampleTest
  tests = functiontests(localfunctions);
end

function test_invariance(testCase)
    ts = TransSyst(3,2);
    ts.add_transition(1,2,1);
    ts.add_transition(2,1,1);
    ts.add_transition(1,3,2);
    ts.add_transition(2,3,2);
    ts.add_transition(3,2,1);
    ts.add_transition(3,3,1);
    ts.add_transition(3,3,2);

    [~, cont1] = ts.pre([1], uint32(1:2), 'exists', 'forall')

    verifyEqual(testCase, cont1(2), uint32(1));

    [~, cont2] = ts.win_always([1, 2], 'exists')

    verifyEqual(testCase, cont2(1), uint32(1));
    verifyEqual(testCase, cont2(2), uint32(1));
end