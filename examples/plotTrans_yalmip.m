%fx = fx1;

mm=1; %mode in trans

syms sx1 sx2

switch mm
    case 1
        fx = [-sx2-0.5*sx1^3-1.5*sx1;...
            -sx2^2+2+sx1];
    case 2
        fx = [-sx2-0.5*sx1^3-1.5*sx1;...
            -sx2+sx1];
    case 3
        fx = [-sx2-0.5*sx1^3-1.5*sx1+2;...
            10+sx1];
    case 4
        fx = [-sx2-0.5*sx1^3-1.5*sx1-1.5;...
            -10+sx1];
end

N_var = 2;
if N_var == 2
    num_pnt_x1 = 20;
    gap = (domain.xmax(1)-domain.xmin(1))/num_pnt_x1;
    [sx1, sx2] = meshgrid(domain.xmin(1):gap:domain.xmax(1),...
    domain.xmin(2):gap:domain.xmax(2));
    dx1 = double(subs(fx(1)));
    dx2 = double(subs(fx(2)));
    mag = sqrt(dx1.^2+dx2.^2);
    dx1u = dx1./mag;
    dx2u = dx2./mag;
    figure; plot(part); hold on; quiver(sx1,sx2,dx1u,dx2u,'k'); axis equal;
end

N_nodes = length(part);
figure; plot(part); hold on;

for i=1:N_nodes
    for j=1:N_nodes
        if (trans(i,j,mm) == 1) && (i~=j)
            p1 = part(i); p2 =part(j);
            xc1 = (p1.xmin+p1.xmax)/2;
            xc2 = (p2.xmin+p2.xmax)/2;
            arrowgit(xc1,xc2,'linewidth',1.5);
        end
    end
end