% Given the state of CTS sys x(t) in Q, map x(t) to [Q]_n
function idx = mapping(xt,Qn,r)
% Input:
%   xt: nx1 vector, state to be mapped
%   Qn: the mesh of discretized state space
% Output:
%   idx: the index of nodes contained by the Br(x(t)) in Qn

if(size(xt,2)~=1)
    error('xt must be a column vector.');
end

if(nargin<2)
    var = load('ts');
    Qn = var.M_X;
    r = var.eta/2;
end

eta = Qn.gridsize;      % grid size of discretized state space
bnd = Qn.bnd;           
n = size(bnd,1);        % # dim
discr_bnd = Qn.discr_bnd;
V = Qn.V; 
% check if xt is out of interested state space
lbnd = (bnd(:,1)<=xt-(r-eta/2));  % lower bnd
ubnd = (xt+(r-eta/2)<=bnd(:,2));  % upper bnd
if((sum(lbnd+ubnd))<2*n)    
    idx=prod(discr_bnd(:,3))+1; % if out of bnd, assign it to sink node
    return;
end

%% find the nodes in Br(xt)
SUB = cell(1,n);       % subscripts of selected discretized nodes
% strN(1:n)=' ';  % param for eval func
for i = 1:n
    pos_node = V{i};%linspace(discr_bnd(i,1),discr_bnd(i,2),discr_bnd(i,3)); % maybe pre-calculation is better
    sub_idx = (abs(pos_node-xt(i))<=(r+1e-15)); % filter the nodes in the ball
                                                % 1e-15 is due to round error in
                                                % matlab
    SUB{i} = find(sub_idx==1);             % find idx           
%     strN(i)=num2str(i);
end

%% map sub to index

X = ndgrid2(SUB);
idx = sub2ind2(discr_bnd(:,3),X);
%  % params of eval func
% str1 (1:(3*n-1))=' ';
% str2 (1:(7*n-1))=' ';
% str1(1:3:(3*n-2))='X';
% str1(2:3:(3*n-1))=strN;
% str1(3:3:(3*n-3))=',';
% str2(1:7:(7*n-6))='S';
% str2(2:7:(7*n-5))='U';
% str2(3:7:(7*n-4))='B';
% str2(4:7:(7*n-3))='{';
% str2(5:7:(7*n-2))=strN;
% str2(6:7:(7*n-1))='}';
% str2(7:7:(7*n-7))=',';
% 
% % generate grid coord X1~Xn
% eval(['[',str1,']=ndgrid2(',str2,');']); 
% 
% % re-use str2
% str2(1:7:(7*n-6))=' ';
% str2(2:7:(7*n-5))='X';
% str2(3:7:(7*n-4))=strN;
% str2(4:7:(7*n-3))='(';
% str2(5:7:(7*n-2))=':';
% str2(6:7:(7*n-1))=')';
% 
% % mapping 
% eval(['idx = sub2ind(discr_bnd(:,3),',str2,');'])


