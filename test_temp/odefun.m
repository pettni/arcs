function dydt=odefun(t,y)
% linear system x = Ax + Bu
persistent A B
if(isempty(A))
%     g = 9.8;
%     h0 = 10;
%     A =[0 1;g/h0 0];
%     B = [0;-1];
    load system
end
dim_x = size(A,2);
dim_u = size(B,2);
dxdt = A*y(1:dim_x) + B*y(dim_x+1:dim_x+dim_u);  % x=y(1:2) u = y(3)
dydt = [dxdt;zeros(dim_u,1)];           % u doesn't change. 