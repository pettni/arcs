% Test calcuation of winning set for specification type
% 'eventually_or_persistence' and the corresponding controller
function tests = exampleTest
  tests = functiontests(localfunctions);
end

function test_win_eventually_or_persistence(testCase)
    ts = TransSyst(7,2);
    ts.add_transition(1,2,1);
    ts.add_transition(2,5,1);
    ts.add_transition(2,6,2);
    ts.add_transition(3,6,1);
    ts.add_transition(4,3,1);
    ts.add_transition(4,2,2);
    ts.add_transition(5,6,1);
    ts.add_transition(6,5,2);
    ts.add_transition(6,7,1);
    ts.add_transition(7,6,2);
    B_list = [5,6,7];
    [W,~,cont]=ts.win_eventually_or_persistence([],{B_list},1);
    verifyEqual(testCase,W,uint32([1,2,3,4,5,6,7]));
end

function test_controller_of_win_eventually(testCase)
    ts = TransSyst(8,2);
    ts.add_transition(1,2,1);
    ts.add_transition(2,1,1);
    ts.add_transition(2,5,1);
    ts.add_transition(3,6,1);
    ts.add_transition(4,3,1);
    ts.add_transition(3,4,1);
    ts.add_transition(4,2,2);
    ts.add_transition(5,6,1);
    ts.add_transition(6,5,2);
    ts.add_transition(6,7,1);
    ts.add_transition(7,6,2);
    ts.add_transition(8,1,1);
    ts.add_transition(8,6,1);
    B_list = [5,6,7];
    [W1,~,cont1]=ts.win_eventually_or_persistence([],{B_list},1);
    verifyFalse(testCase,any(ismember([1,2,3,4],W1)));
    
    ts.add_progress_group([1],[1,2]);
    [W2,~,cont2]=ts.win_eventually_or_persistence([],{B_list},1);
    cont2(8);
    % 8--1-->6 
    cont2(6);
    % 6--2-->5;
    input = cont2(5); % the original one will make error here.
    verifyEqual(testCase,input,uint32(1));
end