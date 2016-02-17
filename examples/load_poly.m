domain = Rec([-2 -1.5; 2 3]);
goal_set = Rec([-1 1.5; -.5 1.9], 1);
unsafe_set = Rec([-2 -1.5; -.5 -1], 2);

sdpvar x1 x2

vars = [x1; x2];

%Moore Greitzer
fx1 = [-x2-0.5*x1^3-1.5*x1;...
    -x2^2+2+x1];

fx2 = [-x2-0.5*x1^3-1.5*x1;...
    -x2+x1];

fx3 = [-x2-0.5*x1^3-1.5*x1+2;...
    10+x1];

fx4 = [-x2-0.5*x1^3-1.5*x1-1.5;...
    -10+x1];

act_set={{fx1},{fx2},{fx3},{fx4}};
%act_set={{fx1},{fx2},{fx3}};

% Build initial partition
part = Partition(domain);
part.add_area(goal_set);
part.add_area(unsafe_set)
part.check();   % sanity check