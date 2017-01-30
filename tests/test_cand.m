%% Main function to generate tests
function tests = exampleTest
    tests = functiontests(localfunctions);
end

function test_forall(testCase)
    ts = TransSyst(4,2);
    ts.add_transition(1,2,1);
    ts.add_transition(2,1,1);

    ts.add_transition(3,1,1);
    ts.add_transition(3,4,1);

    ts.add_transition(1, 3, 2);
    ts.add_transition(2, 3, 2);

    ts.add_transition(3,4,2);
    ts.add_transition(4,3,2);
    ts.add_transition(4,3,1);

    [~, C] = ts.win_always([1, 2, 3], 'exists');
    verifyEqual(testCase, C, [3]);
end 
