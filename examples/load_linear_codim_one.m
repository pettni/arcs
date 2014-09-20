domain = Rec([16 16; 26 26]);
goal_set = Rec([18, 20; 20 22], 1);
unsafe_set = Rec([16 16; 16.2 16.2], 2);

fx1.A = [ -0.002, 0.002; 0 0];
fx1.K = [0;0.1];

fx2.A = [ -0.002, 0.002; 0 0];
fx2.K = [0;-0.1];

fx3.A = [ -0.002, 0.002; 0 0];
fx3.K = [0;0];

fx4.A = [ -0.002, 0; 0 0];
fx4.K = [0.002*16;0];
act_set={fx1,fx2,fx3,fx4};

% Build initial partition
part = Partition(domain);
part.add_area(goal_set);
part.add_area(unsafe_set)
part.check();   % sanity check