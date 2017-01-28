function plot_vf(part, controller)
  % plot_vf(domain, fx, xvar, dvar, dval, color): 
  % plot a 2D SDP vector field over domain
  hold on;
  plot(part, 1, 1);

  xmin = part.domain.xmin(1);
  xmax = part.domain.xmax(1);
  ymin = part.domain.xmin(2);
  ymax = part.domain.xmax(2);
  [X, Y] = meshgrid(linspace(xmin, xmax, 10), linspace(ymin, ymax, 10));

  dx = 0; %(xmax(1) - xmin(1)) / 50;

  color_list = jet(double(part.ts.n_a));

  for a = 1:length(part.act_list)
    fx = part.act_list{a}{1};
    xvar = part.act_list{a}{2};
    dvar = part.act_list{a}{3};
    drec = part.act_list{a}{4};

    fx1 = strrep(strrep(sdisplay(fx(1)), '*', '.*'), '^', '.^');
    fx2 = strrep(strrep(sdisplay(fx(2)), '*', '.*'), '^', '.^');
    fx1 = strrep(fx1, 'xvar(1)', 'x1');
    fx1 = strrep(fx1, 'xvar(2)', 'x2');
    fx1 = strrep(fx1, 'dvar', 'd1');
    fx2 = strrep(fx2, 'xvar(1)', 'x1');
    fx2 = strrep(fx2, 'xvar(2)', 'x2');
    fx2 = strrep(fx2, 'dvar', 'd1');
    eval(strcat('fx1_lambda = @(x1, x2, d1) ', fx1{1}, ';')); 
    eval(strcat('fx2_lambda = @(x1, x2, d1) ', fx2{1}, ';')); 

    for dval = drec.getVertices
      U = fx1_lambda(X+dx*a, Y, dval);
      V = fx2_lambda(X+dx*a, Y, dval);
      quiver(X+dx*a,Y,U,V,'color',color_list(a,:))
    end
  end

  for i=1:length(part)
    if nargin == 2
      a_list = controller(i);
    else
      a_list = 1:part.ts.n_a
    end
    for a=a_list
      for j=part.ts.post(i, a)
        if i ~= j && j <= length(part)-1  % ignore outside state
          [xc1 xc2] = get_best_vector(part.cell_list(i), ...
                                      part.cell_list(j), ...
                                      double(a));
          h = arrow(xc1,xc2,'width',1,'facecolor',color_list(a,:));
        end
      end
    end
  end
end

function [x1, x2] = get_best_vector(rec1, rec2, a)
    % Give endpoints for appropriate arrow
    isect = intersect(rec1, rec2);

    n_dim = isect.getFlatDims();
    t_dim = 3-n_dim;
    xc = getMidpoint(isect);

    xc1 = getMidpoint(rec1);
    xc2 = getMidpoint(rec2);
    dir_n = sign(xc2(n_dim) - xc1(n_dim));
    offs_n = 0.6*min(abs([xc1(n_dim)-xc(n_dim), xc2(n_dim)-xc(n_dim)]));
    
    offs_t = max((isect.xmax - isect.xmin))/2;

    vec_n = zeros(1,2);
    vec_n(n_dim) = dir_n*offs_n;

    vec_t = zeros(1,2);
    vec_t(t_dim) = dir_n*(0.1*offs_t + 0.8*offs_t*(a-1)/5);

    x1 = xc-vec_n+vec_t;
    x2 = xc+vec_n+vec_t;
end
