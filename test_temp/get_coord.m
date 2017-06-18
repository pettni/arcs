function [v1,v2,varargout]=get_coord(X,Mesh)
% input index of nodes and mesh structure, return the coordinates of input
% nodes
nout = max(nargout,1);
x_B = Mesh.ind2sub(X,:); % xt
n = length(Mesh.discr_bnd(:,1));
coord = zeros(n,length(X));

for i = 1:n
    coord(i,:) = Mesh.discr_bnd(i,1)+(x_B(:,i)'-1)*Mesh.gridsize;
end

if(nout>2)
    v1 = coord(1,:)';
    v2 = coord(2,:)';
    for i = 3:nout
        varargout{i-2}=coord(i,:)';
    end
elseif (nout==2)
    v1 = coord(1,:)';
    v2 = coord(2,:)';
else
    v1 = coord(1,:)';
end