%% Main function to generate tests
function tests = exampleTest
	tests = functiontests(localfunctions);
end

%% Test Functions
function test_pre(testCase)

	s = TransSyst(6, 2);
	s.add_transition(3, 1, 1);
	s.add_transition(3, 2, 1);
	s.add_transition(3, 4, 2);

	s.add_transition(4, 1, 1);
	s.add_transition(4, 3, 1);
	s.add_transition(4, 3, 2);
	s.add_transition(5, 4, 1);
	s.add_transition(5, 3, 2);
	s.add_transition(5, 2, 1);
	s.add_transition(5, 1, 2);
	s.add_transition(6, 2, 1);
	s.add_transition(6, 1, 2);
	s.add_transition(6, 2, 2);

	X = [1, 2];
	U = [1, 2];

	verifyEqual(testCase, s.pre(X, U, true, true), uint32([3, 4, 5, 6]))
	verifyEqual(testCase, s.pre(X, U, true, false), uint32([3, 6]))
	verifyEqual(testCase, s.pre(X, U, false, true), uint32([5, 6]))
	verifyEqual(testCase, s.pre(X, U, false, false), uint32([6]))

	s.create_fast()

	verifyEqual(testCase, s.pre(X, U, true, true), uint32([3, 4, 5, 6]))
	verifyEqual(testCase, s.pre(X, U, true, false), uint32([3, 6]))
	verifyEqual(testCase, s.pre(X, U, false, true), uint32([5, 6]))
	verifyEqual(testCase, s.pre(X, U, false, false), uint32([6]))
end

function test_until(testCase)
	s = TransSyst(5, 1);
	s.add_transition(1, 1, 1);
	s.add_transition(1, 2, 1);
	s.add_transition(2, 3, 1);
	s.add_transition(3, 3, 1);
	s.add_transition(4, 3, 1);
	s.add_transition(5, 5, 1);
	s.add_transition(5, 4, 1);

	verifyEqual(testCase, s.win_until(1:5, [3], false), uint32([2, 3, 4]))

	s.add_progress_group([1], [1, 2]);

	verifyEqual(testCase, s.win_until(1:5, [3], false), uint32([1, 2, 3, 4]))

	verifyEqual(testCase, s.win_primal(1:5, [3, 5], {1:5}, 'exists', 'forall'), uint32(1:5));
  verifyEqual(testCase, s.win_dual([], {[]}, [3, 5], 'exists', 'forall'), uint32(1:5));

	s.create_fast()

	verifyEqual(testCase, s.win_until(1:5, [3], false), uint32([1, 2, 3, 4]));

	verifyEqual(testCase, s.win_primal(1:5, [3, 5], {1:5}, 'exists', 'forall'), uint32(1:5));
  verifyEqual(testCase, s.win_dual([], {[]}, [3, 5], 'exists', 'forall'), uint32(1:5));
end

function test_primal1(testCase)
	s = TransSyst(6, 2);

    s.add_transition(1, 2, 1)
    s.add_transition(1, 6, 2)

    s.add_transition(2, 3, 2)
    s.add_transition(2, 6, 1)

    s.add_transition(3, 4, 2)
    s.add_transition(3, 5, 1)

    s.add_transition(4, 1, 1)
    s.add_transition(4, 6, 2)

	s.add_transition(5, 5, 2)
	s.add_transition(5, 6, 1)
	s.add_transition(5, 6, 2)

    % Test []{1 2 3 4 5} and []<>s2 and []<>s1
    safe = [1 2 4 3 5];
    C1 = [3];
    C2 = [2];
    verifyEqual(testCase, ...
    			s.win_primal(safe, 1:6, {C1, C2}, 'exists', 'forall'), ...
                uint32([1 2 3 4]))
    verifyEqual(testCase, ...
    			s.win_primal(safe, 1:6, {C1, C2}, 'forall', 'forall'), ...
                uint32([]))

    s.create_fast()

    verifyEqual(testCase, ...
    			s.win_primal(safe, 1:6, {C1, C2}, 'exists', 'forall'), ...
                uint32([1 2 3 4]))
    verifyEqual(testCase, ...
    			s.win_primal(safe, 1:6, {C1, C2}, 'forall', 'forall'), ...
                uint32([]))

end

function test_primal2(testCase)
	s = TransSyst(5, 2);

	s.add_transition(1, 2, 1);
	s.add_transition(1, 3, 1);
	s.add_transition(1, 4, 1);

	s.add_transition(4, 3, 1);

	s.add_transition(3, 4, 1);
	s.add_transition(3, 2, 1);
	s.add_transition(3, 5, 1);

	s.add_transition(2, 4, 1);

	s.add_transition(5, 3, 1);
	s.add_transition(5, 5, 1);

	c1 = [2];
	c2 = [3, 4];
	tr = 1:5;

	verifyEqual(testCase, ...
				s.win_primal(tr, tr, {c1, c2}, 'exists', 'forall'), ...
				uint32([]))

	s.add_progress_group([1], [3, 4, 5])

	verifyEqual(testCase, ...
				s.win_primal(tr, tr, {c1, c2}, 'exists', 'forall'), ...
				uint32(tr))

    s.create_fast()

	verifyEqual(testCase, ...
				s.win_primal(tr, tr, {c1, c2}, 'exists', 'forall'), ...
				uint32(tr))

end

function test_primal3(testCase)
	s = TransSyst(7, 2);

  s.add_transition(1, 1, 1)
  s.add_transition(1, 1, 2)

  s.add_transition(2, 2, 1)
  s.add_transition(2, 2, 2)

  s.add_transition(3, 1, 1)
  s.add_transition(3, 4, 1)

  s.add_transition(3, 2, 2)
  s.add_transition(3, 3, 2)

  s.add_transition(5, 3, 1)
  s.add_transition(5, 3, 2)
  s.add_transition(5, 5, 2)
  s.add_transition(5, 5, 1)

  s.add_transition(6, 5, 1)
  s.add_transition(6, 6, 1)
  s.add_transition(6, 7, 2)

  s.add_transition(7, 7, 1)
  s.add_transition(7, 6, 1)
  s.add_transition(7, 7, 2)
  s.add_transition(4, 4, 1)
  s.add_transition(4, 4, 2)

  s.add_progress_group([2], [5, 3])

  A = [6, 3, 1];
  B = [3, 2];
  C = [4];

  cC = setdiff(1:7, C);
  cA = setdiff(1:7, A);
  cB = setdiff(1:7, B);
  % win(<>C || <>[]A || <>[]B) = win([]cC && []<>cA && []<>cB )^c
  verifyEqual(testCase, ...
  			s.win_primal(cC, [], {cA, cB}, 'forall', 'forall'), ...
  			uint32([]));

  verifyEqual(testCase, ...
        s.win_dual(C, {A, B}, [], 'exists', 'forall'), ...
        uint32(1:6));

  verifyEqual(testCase, ...
  			s.win_primal(cC, 1:7, {cA, cB}, 'exists', 'forall'), ...
  			uint32(sort(setdiff(1:7, 1:5))));

  verifyEqual(testCase, ...
        s.win_dual(C, {A, B}, [], 'forall', 'forall'), ...
        uint32(1:4));

  s.create_fast()

  verifyEqual(testCase, ...
  			s.win_primal(cC, 1:7, {cA, cB}, 'forall', 'forall'), ...
  			uint32([]));

  verifyEqual(testCase, ...
  			s.win_primal(cC, 1:7, {cA, cB}, 'exists', 'forall'), ...
  			uint32(sort(setdiff(1:7, 1:5))));

end

function test_until_or_always(testCase)

  ts = TransSyst(2,1);
  ts.add_transition(1,1,1);
  ts.add_transition(1,2,1);
  ts.add_transition(2,2,1);

  A = [1];
  P = [2];

  verifyEqual(testCase, ts.win_until_or_always({A}, P, 1), uint32([1 2]));
  verifyEqual(testCase, ts.win_until_or_always({A}, P, 0), uint32([1 2]));

  ts = TransSyst(7,1);
  ts.add_transition(1,3,1);
  ts.add_transition(2,3,1);
  ts.add_transition(3,4,1);
  ts.add_transition(4,5,1);
  ts.add_transition(5,5,1);
  ts.add_transition(4,6,1);
  ts.add_transition(6,7,1);

  A1 = [1, 3, 4, 6];
  A2 = [2, 3, 4, 5];
  P = [7];

  verifyEqual(testCase, ts.win_until_or_always({A1, A2}, P, 1), uint32([3,4,5,6,7]));

  ts = TransSyst(7,2);
  ts.add_transition(1,3,1);
  ts.add_transition(2,3,1);
  ts.add_transition(3,4,1);
  ts.add_transition(4,5,2);
  ts.add_transition(5,5,1);
  ts.add_transition(4,6,1);
  ts.add_transition(6,7,1);
  ts.add_transition(7,7,1);

  A1 = [1, 3, 4, 6, 7];
  A2 = [2, 3, 4, 5];
  P = [];

  verifyEqual(testCase, ts.win_until_or_always({A1, A2}, P, 1), uint32([1,2,3,4,5,6,7]));

  [~, ~, cont] = ts.win_until_or_always({A1, A2}, P, 1);
  state = uint32([3]);
  for i=1:20
    act = cont(state(end));
    state(end+1) = ts.post(state(end), act(randi(length(act))));
  end
  verifyEqual(testCase, or(all(ismember(state, A1)), all(ismember(state, A2))), true);
end