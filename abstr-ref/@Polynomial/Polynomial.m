classdef Polynomial<handle
  % Create a Polynomial
  % 
  % - Polynomial(expr, vars) create a Polynomial from sdpvar expression
  %                          'expr' in variables 'var'.
  %
  % - Polynomial(n) create the zero Polynomial in n variables
  %
  % - Polynomial(n, v) create the Polynomial in n variables with grlex 
  %                    coefficient vector 'v'
  %
  % - Polynomial(coef, mons) create a Polynomial with coefficients 'coef' (1 x k)
  %                          for the monomials 'mons' (n x k), for k > 1
  %
  % Remark: this class is experimental and in need of improvement for general-
  %         purpose use

  properties (SetAccess=protected)
    coef;
    mons;
  end

  methods
    function p = Polynomial(arg1, arg2)
      if nargin > 1 && isa(arg2, 'sdpvar')
        % sdpvar case
        if length(arg1) > 1
          p = [Polynomial(arg1(1), arg2)];
          for i=2:length(arg1)
            p(end+1) = Polynomial(arg1(i), arg2);
          end
          return
        end
        [coef, mono_text] = coefficients(arg1);
        p.coef = coef;
        p.mons = zeros(length(arg2), length(coef));
        for i=1:length(coef)
          p.mons(:,i) = degree(mono_text(i), arg2);
        end
      elseif length(arg1) == 1
        % grlex vector
        p.coef = zeros(0,0);
        p.mons = zeros(arg1, 0);
        p.coef = [];
        p.mons = zeros(arg1, 0);
        if nargin == 2
          for i = 1:length(arg2)
            p.coef(end+1) = arg2(i);
            p.mons(:, end+1) = mono_unrank_grlex(arg1, i);
          end
        end
      else
        % monomial/coefficient vector
        p.coef = arg1;
        p.mons = arg2;
      end
    end

    function reduce(p)
      % Find and sum doubles
      [~, iproj, ilift] = unique(p.mons', 'rows');
      p.mons = p.mons(:, iproj);
      newcoef = zeros(1, length(iproj));
      for i=1:length(iproj)
        newcoef(i) = sum(p.coef(ilift == i));
      end
      p.coef = newcoef;

      % Remove zero coefficients
      p.mons = p.mons(:, p.coef ~= 0);
      p.coef = p.coef(p.coef ~= 0);
    end

    function clean(p, eps)
      % Remove terms smaller than 'eps'
      if nargin < 2
        eps = 1e-6;
      end
      p.mons = p.mons(:, abs(p.coef) > eps);
      p.coef = p.coef(:, abs(p.coef) > eps);
    end

    function ret = dim(p)
      % return the number of variables
      ret = size(p.mons, 1);
    end

    function ret = deg(p)
      % return the degree
      if length(p) > 1
        ret = 0;
        for i = 1:length(p)
          ret = max(ret, p(i).deg);
        end
        return
      end
      ret = max(sum(p.mons, 1));
    end

    function ret = mon_vec(p, deg)
      % Return the coefficient vector in the grlex monomial basis
      if nargin == 1
        ret = zeros(count_monomials_leq(p.dim(), p.deg()), 1);
      else
        if deg < p.deg
          error('degree to low')
        end
        ret = zeros(count_monomials_leq(p.dim(), deg), 1);
      end
      for i=1:length(p.coef)
        ret(mono_rank_grlex(p.dim, p.mons(:,i))) = p.coef(i);
      end
    end

    function ret = evaluate(p, point)
      % Evaluate the polynomial at 'points' (n \times k)
      if size(point, 2) == 1
        ret = sum(p.coef .* prod(point.^p.mons, 1));
      else
        ret = zeros(1, size(point, 2));
        for i=1:size(point,2)
          ret(1,i) = p.evaluate(point(:,i));
        end
      end
    end

    function disp(p)
      if length(p) > 1
        disp('Vector of Polynomial')
        return
      end
      disp(['Polynomial in ', num2str(p.dim), ' variables:'])
      str = cell(1, length(p.coef));
      for i = 1:length(p.coef)
        str{i} = strcat(num2str(p.coef(i)), ' * x^(', num2str(p.mons(:,i)'), ')');
      end
      disp(strjoin(str, ' + '));
    end

    function ret = plus(p1, p2)
      if p1.dim ~= p2.dim
        error('dimension mismatch')
      end
      ret = Polynomial(p1.dim);
      ret.coef = [p1.coef p2.coef];
      ret.mons = [p1.mons p2.mons];
      ret.reduce
    end

    function ret = uminus(p)
      ret = Polynomial(-p.coef, p.mons);
    end

    function ret = minus(p1, p2)
      ret = plus(p1, -p2);
    end

    function ret = mtimes(p1, p2)
      ret = Polynomial(p1.dim);
      if isa(p2, 'Polynomial')
        for i = 1:length(p1.coef)
          for j = 1:length(p2.coef)
            ret.coef(end+1) = p1.coef(i) * p2.coef(j);
            ret.mons(:, end+1) = p1.mons(:,i) + p2.mons(:,j);
          end
        end
      else
        ret.coef = p2 * p1.coef;
        ret.mons = p1.mons;
      end
      ret.reduce;
    end

    function q = shift(p, xs)
      % Compute the shifted polynomial q(x) = p(x+xs)
      q = Polynomial(p.dim);
      xs = reshape(xs, length(xs), 1);
      for term = 1:length(p.coef)
        q_i = Polynomial(p.dim, [1]);
        for j = 1:p.dim
          alpha_ij = p.mons(j, term);
          q_ij = Polynomial(p.dim);
          q_ij.coef = arrayfun(@(k) nchoosek(alpha_ij,k), 0:alpha_ij).*(xs(j).^(alpha_ij:-1:0));
          q_ij.mons = zeros(p.dim, alpha_ij+1);
          q_ij.mons(j, :) = 0:alpha_ij;
          q_i = q_i * q_ij;
        end
        q = q + q_i * p.coef(term);
      end
      q.reduce;
    end

    function q = mrdivide(p, d)
      % Divide a polynomial by a double d
      q = Polynomial(p.dim);
      q.coef = p.coef/d;
      q.mons = p.mons;
    end

    function q = scale(p, xs)
      % Compute the scaled polynomial q(x) = p(xs.*x)
      q = Polynomial(p.dim);
      q.coef = p.coef;
      q.mons = p.mons;
      xs = reshape(xs, length(xs), 1);
      for term=1:length(p.coef)
        q.coef(term) = q.coef(term) * prod(xs.^q.mons(:,term));
      end
    end

    function res = eq(p1, p2)
      if length(p1.coef) ~= length(p2.coef)
        res = false;
        return
      end
      p1.reduce;  % sorts
      p2.reduce;  % sorts
      res = norm(p1.coef-p2.coef) < 1e-5 && all(all(p1.mons == p2.mons));
    end

    function res = isequal(p1, p2)
      res = p1 == p2;
    end

  end

end