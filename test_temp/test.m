clear all;clc;
addpath(genpath('../abstr-ref'));
load test.mat
%% Test 'win_eventually_or_persistence'

disp("Test 'win_eventually_or_persistence'")
pause(0.1);
ts.create_fast();
[W, C, cont]=ts.win_eventually_or_persistence([],{B_list'},1);

% Visualization
visual
%% Test 'Controller'

disp("Test 'Controller'")

x0 = [-2.43;3.5]; % initial condition
idx_x0 = mapping(x0,M_X,eta/2);
if(~ismember(idx_x0,W))
    error('x0 is beyond winning set.');
end
[x1,x2] =get_coord(idx_x0,M_X);
plot(x1,x2,'xb','markersize',10)

idx_x = idx_x0;
[x1,x2] = get_coord(idx_x,M_X);
xt = [x0];
idx_u = 0;

X_list = [idx_x];

t_span = 60;
for i = 1:t_span
    % visual (on the grid)
    % get the options of input 
%     disp(idx_x)
    u_option = cont.get_input(idx_x);
    
    idx_u = u_option(1);
     
    u0 = get_coord(idx_u,M_U);    % get the coordinate of input
    y0 = [xt;u0];    % States for numerical integration
    
    yt = ode45(@odefun,[0,tau],y0);
    
    xt = yt.y(1:2,end); % destination in one step
    idx_x  = mapping(xt,M_X,eta/2);
    
    X_list = [X_list;idx_x];
    figure(2);
    plot((i-1)*tau+yt.x,yt.y(1,:),'r-'); % position red
    plot((i-1)*tau+yt.x,yt.y(2,:),'b-'); % velocity blue
    hold on;
    drawnow;
end

legend('Vel of Mass Center','Pos of Mass Center');
xlabel('t');

figure(1);

[x1,x2] = get_coord(X_list,M_X);
for i=1:length(x1)-1
    arrow('Start',[x1(i),x2(i)],'Stop',[x1(i+1),x2(i+1)],'Length',10,'TipAngle',5)
    pause(0.05)
end

disp('Done.');