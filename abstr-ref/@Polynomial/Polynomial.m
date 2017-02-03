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

  properties (SetAccess=protected)
    coef;
    mons;
  end

  methods
    function p = Polynomial(arg1, arg2)
      if isa(arg1, 'sdpvar')
        % sdpvar case
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
        for i = 1:length(arg2)
          p.coef(end+1) = arg2(i);
          p.mons(:, end+1) = mono_unrank_grlex(arg1, i);
        end
      else
        % monomial/coefficient vector
        p.coef = arg1;
        p.mons = arg2;
      end
      p.reduce;
    end

    function reduce(p)
      % Remove zero coefficients
      p.mons = p.mons(:, p.coef ~= 0);
      p.coef = p.coef(p.coef ~= 0);
    end

    function add_term(p, coef, mon)
      % Add a term to the monomial
      p.coef(end+1) = coef;
      p.mons(:, end+1) = mon;
    end

    function ret = dim(p)
      % return the number of variables
      ret = size(p.mons, 1);
    end

    function ret = deg(p)
      % return the degree
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

  end

end