classdef PolyLinTrans<handle
  properties (SetAccess=protected)
    n0 = [];   % rhs number of vars  
    n1 = [];   % lhs number of vars
    d0 = [];   % rhs degree
    d1 = [];   % lhs degree
    data = {};   % row_idx, col_idx, coeff
  end

  methods
    function plt = PolyLinTrans(n0, n1, d0, d1)
      % Linear transformation of a polynomial in the monomial basis.
      plt.n0 = n0;
      plt.n1 = n1;
      plt.d0 = d0;
      plt.d1 = d1;
      plt.data = sparse(count_monomials_leq(n0, d0), ...
                        count_monomials_leq(n1, d1));
    end

    function varargout = subsref(plt, S)
      switch S(1).type
        % call builtin method to access class properties
      case '.' 
        [varargout{1:nargout}] = builtin('subsref', plt, S);
        % delegate to cell_list indexing
      case '()'
        try
          varargout{1} = plt.data(S.subs{1}, S.subs{2});
        catch
          varargout{1} = 0;
        end
      end 
    end

    function plt = plus(plt1, plt2)
      if plt1.n0 ~= plt2.n0 || plt1.n1 ~= plt2.n1
        error('dimension mismatch in plus')
      end
      plt = PolyLinTrans(plt1.n0, plt1.n1, max(plt1.d0, plt2.d0), ...
                                           max(plt1.d1, plt2.d1));
      
      [col1, row1, val1] = find(plt1.data);
      [col2, row2, val2] = find(plt2.data);
      
      ncol = count_monomials_leq(plt.n1, plt.d1);
      nrow = count_monomials_leq(plt.n0, plt.d0);

      plt.data = sparse(col1, row1, val1, nrow, ncol) + ...
                 sparse(col2, row2, val2, nrow, ncol);

    end

    function ret = mtimes(T, p)
      % Compute T p for a PolyLinTrans T and polynomial p
      % or for two PolyLinTrans
      if isa(p, 'Polynomial')
        if T.n0 ~= p.dim
          error('dimension mismatch')
        end
        if T.d0 < p.deg
          error('degree is too high')
        end
        ret_monv = T.as_vector_trans * p.mon_vec;
        ret = Polynomial(T.n1, ret_monv);

      elseif isa(p, 'PolyLinTrans')
        if p.n1 ~= T.n0
          error('dimension mismatch');
        end
        if p.d1 > T.d0
          error('degree low')
        end

        ret = PolyLinTrans(p.n0, T.n1, p.d0, T.d1);
        ret.data = p.data * T.data(1:count_monomials_leq(p.n1, p.d1), :);
      end
    end

    function ret = as_vector_trans(plt)
      % Return a representation A of the transformation from a vector
      % representing a symmetric matrix.
      %
      % That is, if v is a monomial coefficient vector for p = v' * m(x), then 
      %    A v 
      % is the monomial coefficient vector of the transformed polynomial 
      % T p = (Av)' * m(x)
      ret = plt.data';
    end

    function ret = as_matrix_trans(plt)
      % Return a representation A of the transformation from a vector
      % representing a symmetric matrix.
      %
      % That is, if v = vec(S) represents the symmetric matrix S, then
      %    A v 
      % is the monomial coefficient vector of the transformed polynomial,
      % i.e. m(x,d/2)' * S * m(x,d/2) = (A v)' * m(x,d), where d is degree

      half_deg = ceil(plt.d0/2);
      num_mon = count_monomials_leq(plt.n0, half_deg);
      L = num_mon * (num_mon + 1) / 2;  % length of vector
      n = (sqrt(1 + 8 * L) - 1)/2;  % side of matrix

      j_mat = 1;
      i_mat = 1;
      grlex_i = zeros(plt.n0, 1);
      grlex_j = zeros(plt.n0, 1);

      ret = sparse(count_monomials_leq(plt.n1, plt.d1), L);

      for k_vec = 1:L
        % Step through matrix
        col_idx = mono_rank_grlex(plt.n0, grlex_i + grlex_j);
        try
          [~, row_idx, val] = find(plt.data(col_idx, :));
          ret(row_idx, k_vec) = val * (2 - (i_mat == j_mat));
        end
        
        if j_mat == n
          i_mat = i_mat+1;
          j_mat = i_mat;
          grlex_i = mono_next_grlex(plt.n0, grlex_i);
          grlex_j = grlex_i;
        else
          grlex_j = mono_next_grlex(plt.n0, grlex_j);
          j_mat = j_mat+1;
        end
      end
    end

    function disp(plt, all) 
      if nargin < 2
        all = false;
      end
      disp(['Polynomial transformation P(', num2str(plt.n0), ',', num2str(plt.d0), ') --> P(', num2str(plt.n1), ',', num2str(plt.d1), ')'])
      [row,col,val] = find(plt.data);
      if all
        for i = 1:length(row)
          from = mono_unrank_grlex(plt.n0, row(i));
          to = mono_unrank_grlex(plt.n1, col(i));
          disp(['x^(', num2str(from'), ') |--> ', num2str(val(i)), ' * x^(', num2str(to'), ')']);
        end
      end
    end
  end

  methods (Static)

    function plt = eye(n0, n1, d0, d1)
      % Identity transformation
      if n1 < n0
        error('eye requires n1 >= n0')
      end
      plt = PolyLinTrans(n0, n1, d0, d1);
      if n0 == n1
        for i = 1:count_monomials_leq(n0, d0)
          plt.data(i, i) = 1.;
        end
      else
        for i = 1:count_monomials_leq(n0, d0)
          idx_1 = mono_unrank_grlex(n0, i);
          j = mono_rank_grlex(n1, [idx_1; zeros(n1-n0, 1)]);
          plt.data(i, j) = 1.;
        end
      end
    end

    function plt = diff(n, d, xi)
      % Transform p |-> d_xi p
      % Lowers degree by one
      plt = PolyLinTrans(n, n, d, d-1);
      for i = 1:count_monomials_leq(n, d)
        grlex_i = mono_unrank_grlex(n, i);
        k = grlex_i(xi);
        if k > 0
          grlex_i(xi) = k - 1;
          plt.data(i, mono_rank_grlex(n, grlex_i)) = k;
        end
      end
    end

    function plt = elvar(n, d, xi, val)
      % Transform p(x(i-1), xi, x(i+1)) |-> p(x(i-1), val, x(i+1))
      % Decreases number of variables by 1
      if length(xi) ~= length(val)
        error('xi and val must be same length')
      end
      plt = PolyLinTrans(n, n-length(xi), d, d);
      val = reshape(val, length(val), 1);
      for i = 1:count_monomials_leq(n, d)
        grlex_i = mono_unrank_grlex(n, i);
        coeff = prod(val.^grlex_i(xi));
        grlex_i(xi) = [];
        j = mono_rank_grlex(n-length(xi), grlex_i);
        plt.data(i, j) = coeff;
      end
    end

    function plt = mul_pol(n, d0, d1, q)
      % mul_pol(n, d0, d1, q): Transformation p |-> q * p
      % d0: initial degree
      % d1: final degree: must be >= d0 + deg(q)
      % q : Polynomial

      if ~isa(q, 'Polynomial') && q.dim ~= n
        error('variable mismatch')
      end
      if d0 + q.deg > d1
        error('degree too low')
      end

      q.reduce();  % get rid of redundant terms

      plt = PolyLinTrans(n, n, d0, d1);
      for j = 1:size(q.mons, 2)  % for each term
        grlex_i = zeros(n, 1);
        i = 1;
        while sum(grlex_i) <= d0   % for each monomial
          plt.data(i, mono_rank_grlex(n, grlex_i + q.mons(:, j))) = q.coef(j);
          grlex_i = mono_next_grlex(n, grlex_i);
          i = i+1;
        end
      end

    end
  end

end