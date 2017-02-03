function tests = exampleTest
  tests = functiontests(localfunctions);
end

function test_diff(testCase)
	L = PolyLinTrans.diff(2,2,1);

	verifyEqual(testCase, L(3,1), sparse(1));
	verifyEqual(testCase, L(5,2), sparse(1));
	verifyEqual(testCase, L(6,3), sparse(2));
end

function test_elvar(testCase)
  L = PolyLinTrans.elvar(2, 2, 1, 0);

  verifyEqual(testCase, L(1,1), sparse(1));
  verifyEqual(testCase, L(2,2), sparse(1));
  verifyEqual(testCase, L(4,3), sparse(1));

  L = PolyLinTrans.elvar(3,3,[1,2],[0,3]);

  verifyEqual(testCase, L(1,1), sparse(1));
  verifyEqual(testCase, L(2,2), sparse(1));
  verifyEqual(testCase, L(3,1), sparse(3));
  verifyEqual(testCase, L(4,1), sparse(0));
  verifyEqual(testCase, L(9,1), sparse(0));
  verifyEqual(testCase, L(3,1), sparse(3));
  verifyEqual(testCase, L(13,2), sparse(9));
end

function test_mulpol(testCase)
  grlex = [2  0 1;
           0  2 0];
  coefs = [1 -1 3];
  L = PolyLinTrans.mul_pol(2, 3, Polynomial(coefs, grlex));

  verifyEqual(testCase, L.d0, 3);
  verifyEqual(testCase, L.d1, 5);

  verifyEqual(testCase, L(mono_rank_grlex(2, [0,0]'), ...
                          mono_rank_grlex(2, [1,0]')), sparse(3));
end

function test_as_matrix_trans(testCase)
  L = PolyLinTrans.eye(2, 2, 4, 4);

  vecA = 1:21;
  A = [1, 2, 3,   4,  5,  6;
        2, 7, 8,   9, 10, 11;
        3, 8, 12, 13, 14, 15;
        4, 9, 13, 16, 17, 18;
        5, 10,14, 17, 19, 20;
        6, 11,15, 18, 20, 21];

  % Up to degree 4 for p' mon1
  mon1c = [0 0 1 0 1 2 0 1 2 3 0 1 2 3 4;
           0 1 0 2 1 0 3 2 1 0 4 3 2 1 0];

  % Up to degree 2 for mon2' A mon2
  mon2c = [0 0 1 0 1 2;
           0 1 0 2 1 0];

  % Should have mon2' A mon2 = (L vecA)' mon1
  for i=1:10
    x = randn(1,1);
    y = randn(1,1);
    mon1 = prod([x;y].^mon1c, 1);
    mon2 = prod([x;y].^mon2c, 1);
    AL = L.as_matrix_trans();
    verifyEqual(testCase, mon2 * A * mon2', ...
                mon1 * AL * vecA', ...
                'reltol', 1e-10);
  end
end

function test_as_vector_trans(testCase)
  grlex = [2  0 1;
           0  2 0];
  coefs = [1 -1 3];
  L = PolyLinTrans.mul_pol(2, 3, Polynomial(coefs, grlex));
  AL = L.as_vector_trans();
  for i = 1:size(AL,1)
    for j = 1:size(AL,2)
      verifyEqual(testCase, AL(i,j), sparse(L(j,i)));
    end
  end
end

function test_ij_k(testCase)
  verifyEqual(testCase, symmat_ij_k(1,1,21), 1);
  verifyEqual(testCase, symmat_ij_k(1,2,21), 2);
  verifyEqual(testCase, symmat_ij_k(2,1,21), 2);
  verifyEqual(testCase, symmat_ij_k(6,6,21), 21);
  verifyEqual(testCase, symmat_ij_k(1,6,21), 6);
  verifyEqual(testCase, symmat_ij_k(2,2,21), 7);
end

function test_k_ij(testCase)
  [i,j] = symmat_k_ij(1, 21);
  verifyEqual(testCase, i, 1);
  verifyEqual(testCase, j, 1);

  [i,j] = symmat_k_ij(2, 21);
  verifyEqual(testCase, i, 1);
  verifyEqual(testCase, j, 2);

  [i,j] = symmat_k_ij(2, 21);
  verifyEqual(testCase, i, 1);
  verifyEqual(testCase, j, 2);

  [i,j] = symmat_k_ij(21, 21);
  verifyEqual(testCase, i, 6);
  verifyEqual(testCase, j, 6);

  [i,j] = symmat_k_ij(6, 21);
  verifyEqual(testCase, i, 1);
  verifyEqual(testCase, j, 6);

  [i,j] = symmat_k_ij(7, 21);
  verifyEqual(testCase, i, 2);
  verifyEqual(testCase, j, 2);
end