domain = Rec([-150, -150; 150, 150]);
goal_set = Rec([-140, 50; -50 140], 1);
unsafe_set = Rec([-150 140; 150 150], 2);

output_impedence = 1;
load1 = 10;
load2 = 100000;
V_ref = 110;

Kp1 = 2;
Ki1 = 1;

c1 = load1/(output_impedence+load1);
c2 = load2/(output_impedence+load2);

%load 1
fx1.A = [-Kp1*c1 Kp1-Ki1; c1 -1];
fx1.K = [-Kp1^2*c1+Ki1; Kp1]*V_ref;

%load 2
fx2.A = [-Kp1*c2 Kp1-Ki1; c2 -1];
fx2.K = [-Kp1^2*c2+Ki1; Kp1]*V_ref;

act_set={fx1,fx2};

% Build initial partition
part = Partition(domain);
part.add_area(goal_set);
part.add_area(unsafe_set)
part.check();   % sanity check
