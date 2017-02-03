classdef Rec<handle

  properties (SetAccess=protected)
    dim;    % Dimension of rec, scalar
    xmin;   % Row vector of lower bounds
    xmax;   % Row vector of upper bounds
    ap;     % List of atomic propositions in cell
  end

  methods
    function r = Rec(X, ap)
      % Rec: Create a hyperrectangle
      % 
      % SYNTAX
      % ------
      %
      %   rec = Rec(X);
      %   rec = Rec(X,ap);
      %
      % INPUT
      % -----
      %   
      %   X   n x 2 matrix defining min and max along each dimension.
      %   ap  cell of APs of of the hyperrectangle. Default: {}
      %
      % PROPERTIES
      % ------
      %
      %   dim     dimension
      %   xmin    array (1 x dim) of lower bounds for each dimension
      %   xmax    array (1 x dim) of upper bounds for each dimension
      % 
      % SEE ALSO
      % -------
      % 
      % add_ap, remove_ap, isFullDim, isEmptySet, getFlatDims, getFullDims, getMidpoint, 
      % volume, split, isInside, contains, contains_strictly, getVertices, getVertexI,
      % getFacetI, isNeighbor, mldivide, plot, projection, toPoly
      if nargin<2
        ap = {};
      end

      if size(X,1) ~= 2
        X = X';
      end
      if size(X,1) ~= 2
        error('Wrong dimensions of X')
      end
      r.dim = size(X,2);
      r.xmin = X(1,:);
      r.xmax = X(2,:);
      r.ap = ap;
    end

    function add_ap(rec, ap)
      % Add an atomic proposition
      rec.ap = union(rec.ap, ap);
    end

    function remove_ap(rec, ap)
      % Remove an atomic proposition
      rec.ap = setdiff(rec.ap, ap);
    end

    function ret = isFullDim(rec)
      % Returns true if Rec is fully dimensional
      ret = all(rec.xmin<rec.xmax);
    end

    function ind = getFlatDims(rec)
      % Returns indices of "flat" dimensionssuch that x_{i,min} = x_{i,max}
      if rec.isEmptySet
        ind = [];
        return;
      end
      ind = find(rec.xmin == rec.xmax);
    end

    function ind = getFullDims(rec)
      % Returns indices of "full" dimensions i such that x_{i,min} < x_{i,max}
      if rec.isEmptySet
        ind = [];
        return;
      end
      ind = find(rec.xmin < rec.xmax);
    end

    function ret = isEmptySet(rec)
      % Returns true if Rec is the empty set
      if length(rec) > 1
        ret = false(1,length(rec));
        for i=1:length(rec)
          ret(i) = isEmptySet(rec(i));
        end
        return
      end 
      ret = any(rec.xmin>rec.xmax);
    end

    function mid = getMidpoint(rec)
      % Returns the middle point of a Rec 
      mid = (rec.xmin + rec.xmax)/2;
    end

    function vol = volume(rec)
      % Computes volume
      if isempty(rec)
        vol = 0;
        return
      end

      if length(rec)>1
        vol = zeros(1,length(rec));
        for i = 1:length(rec)
          vol(i) = volume(rec(i));
        end
        return;
      end

      vol = prod(rec.xmax-rec.xmin);
    end

    function [part1, part2] = split(rec, dim)
      % Split the Rec along dimension dim
      % if not given, choose maximal
      if nargin<2
        [~, dim] = max(rec.xmax - rec.xmin);
      end
      xmin = rec.xmin;
      xmax = rec.xmax;
      midval = (xmax(dim) + xmin(dim))/2;
      xminspl = xmin; xminspl(dim) = midval;
      xmaxspl = xmax; xmaxspl(dim) = midval;
      part1 = Rec([xmin; xmaxspl], rec.ap);
      part2 = Rec([xminspl; xmax], rec.ap);
    end

    function ret = eq(rec1, rec2)
      % Overloaded == operator
      if length(rec1)>1 || length(rec2)>1
        ret1 = mldivide(rec1, rec2);
        ret2 = mldivide(rec2, rec1);
        ret = isempty(ret1) && isempty(ret2);
        return;
      end
      ret = all(rec1.xmin == rec2.xmin) && all(rec1.xmax == rec2.xmax);
    end

    function ret = ne(rec1, rec2)
      % Overloaded ~= operator
      ret = ~eq(rec1,rec2);
    end

    function ret = isInside(rec, point)
      % Return true if point is inside rec
      if length(rec) > 1
        ret = false;
        for i=1:length(rec)
          if isInside(rec(i), point)
            ret = true;
          end
        end
        return
      end

      if length(point) ~= rec.dim
        error('isInside: dimension mismatch')
      end

      point = reshape(point,1,length(point));
      ret = all(rec.xmin<=point) && all(point<=rec.xmax);
    end

    function ret = contains(rec1, rec2)
      % returns true if rec2 is contained in rec1
      ret = all(rec2.xmax<=rec1.xmax) && all(rec2.xmin>=rec1.xmin);
    end

    function ret = contains_strictly(rec1, rec2)
      % Return true if rec2 does not touch the boundaries of rec1
      ret = all(rec1.xmax < rec2.xmax) && all(rec1.xmin>rec2.xmin);
    end

    function p = bounding_polynomial(rec, i)
      % Return the Polynomial p = (x_i^+ - x)*(x - x_i^-)
      dm = rec.xmin(i);
      dp = rec.xmax(i);
      mons = zeros(rec.dim, 3);
      mons(i, :) = [0 1 2];
      p = Polynomial([-dm*dp dm+dp -1], mons);
    end

    function r = elvar(rec, i)
      % return a lower-dimensional Rec where the i'th dimension 
      % has been projected away
      xm = rec.xmin(1:rec.dim ~= i);
      xp = rec.xmax(1:rec.dim ~= i);
      r = Rec([xm; xp]);
      r.ap = rec.ap;
    end

    function r = mtimes(rec1, rec2)
      % return higher-dimensional Rec rec1 x rec2
      r = Rec([rec1.xmin rec2.xmin; rec1.xmax rec2.xmax]);
    end

  end
end