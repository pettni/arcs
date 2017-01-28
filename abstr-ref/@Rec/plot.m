function handle = plot(rec, color, alpha, multiple)
  % Rec/PLOT: Plot a rectangle or an array of rectangles.
  % 
  % SYNTAX
  % ------
  %
  % h = plot(rec);
  % h = plot(rec,color);
  % h = plot(rec,color,alpha);
  % h = plot(rec,color,alpha,mult);
  % 
  % INPUT
  % -----
  % 
  % color   color to use in the form of a 1x3 RGB array. Default [1 0 0] (red).
  %     Class: double
  %   alpha   transparency level between 0 and 1 (useful in 3D). Default 1.
  %   Class: double
  % mult  plot rectangles in the array in different colors. Default 1.
  %   Class: bool
  %
  % OUTPUT
  % ------
  %
  %   h   handle of graphics object. In the case of an array, 
  %   the handle of the last rectangle is returned (useful for labels).
  %   Class: graphics handle
  if nargin<2
    color = [1 0 0];
  end
  if nargin<3
    alpha = 1;
  end
  if nargin<4
    multiple = 1;
  end
  if length(rec)>1
    % We have an array
    hh = ishold; if (~hh) clf; end
    hold on
    if multiple
      color = autumn(length(rec));
    end
    for i = 1:length(rec)
      handle = plot(rec(i), color(min(i, size(color,1)),:), alpha);
    end
    if (~hh) hold off; end
    return;
  end 

  if length(rec) == 0
    disp('Warning: tried to plot nonexistent cell')
    % do nothing
    return;
  end

  if rec.dim>3
    error('Cant plot in dimension larger than 3')
  end

  if rec.isEmptySet
    return;
  end

  hh = ishold;
  if ~hh
    clf
  end
  hold on;

  if rec.dim==1
    h = line([rec.xmin rec.xmax], [0 0]);
    set(h, 'Color', color);
  elseif rec.dim==2
    if rec.isFullDim
      vert = rec.getVertices;
      h = fill(vert(:,1), vert(:,2), 'r');
      set(h, 'FaceColor', color);
      set(h, 'FaceAlpha', alpha);
    else
      h = line(vert(:,1), vert(:,2), 'r');
      set(h, 'Color', color);
    end
  elseif rec.dim==3
    n_fulldim = length(rec.getFullDims);
    if n_fulldim == 3
      projx = projection(rec, [2 3]);
      vertx = projx.getVertices;
      h = fill3(rec.xmin(1)*ones(4,1), vertx(:,1), vertx(:,2), 'r'); set(h, 'FaceColor', color); set(h, 'FaceAlpha', alpha)
      h = fill3(rec.xmax(1)*ones(4,1), vertx(:,1), vertx(:,2), 'r'); set(h, 'FaceColor', color); set(h, 'FaceAlpha', alpha)
      projy = projection(rec, [1 3]);
      verty = projy.getVertices;
      h = fill3(verty(:,1), rec.xmin(2)*ones(4,1), verty(:,2), 'r'); set(h, 'FaceColor', color); set(h, 'FaceAlpha', alpha)
      h = fill3(verty(:,1), rec.xmax(2)*ones(4,1), verty(:,2), 'r'); set(h, 'FaceColor', color); set(h, 'FaceAlpha', alpha)
      projz = projection(rec, [1 2]);
      vertz = projz.getVertices;
      h = fill3(vertz(:,1), vertz(:,2), rec.xmin(3)*ones(4,1), 'r'); set(h, 'FaceColor', color); set(h, 'FaceAlpha', alpha)
      h = fill3(vertz(:,1), vertz(:,2), rec.xmax(3)*ones(4,1), 'r'); set(h, 'FaceColor', color); set(h, 'FaceAlpha', alpha)
    elseif n_fulldim == 2
      ind_fulldims = rec.getFullDims;
      ind_flatdim = rec.getFlatDims;
      flatdimval = rec.xmin(ind_flatdim);
      proj = projection(rec, ind_fulldims);
      vert = getVertices(proj);
      if ind_flatdim == 1
        h = fill3(flatdimval*ones(4,1), vert(:,1), vert(:,2), 'r');
      elseif ind_flatdim == 2
        h = fill3(vert(:,1), flatdimval*ones(4,1), vert(:,2), 'r');
      else
        h = fill3(vert(:,1), vert(:,2), flatdimval*ones(4,1), 'r');
      end
      set(h, 'FaceColor', color)
      set(h, 'FaceAlpha', alpha)
    else
      vert = rec.getVertices;
      h = line(vert(:,1), vert(:,2), vert(:,3));
      set(h, 'Color', color)
    end
  end
  if ~hh
    hold off
  end
  handle = h;
end