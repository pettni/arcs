clear all
yalmip clear

global ops
ops = sdpsettings('solver', 'mosek', 'cachesolvers', 1, 'verbose', 0);

domain = Rec([-2 -1.5; 2 3]);
goal_set = Rec([-1 0.5; -0.2 1.8], {'goal'});
unsafe_set = Rec([-2 -1.5; -.5 -1], {'unsafe'});
d_rec = Rec([-5; 5]);

domain = Rec([-2 -1.5; 2 3]);
sdpvar x1 x2 d1

%Moore Greitzer
fx = {};
fx{1} = [-x2-0.5*x1^3-1.5*x1;...
        -x2^2+2+x1+d1];
fx{2} = [-x2-0.5*x1^3-1.5*x1;...
        -x2+x1+d1];
fx{3} = [-x2-0.5*x1^3-1.5*x1+2;...
        10+x1+d1];
fx{4} = [-x2-0.5*x1^3-1.5*x1-1.5;...
        -10+x1+d1];

part = Partition(domain);
part.add_area(goal_set);
part.add_area(unsafe_set)
part.create_ts();

part.add_mode({fx{1}, [x1; x2], [d1], d_rec});
part.add_mode({fx{2}, [x1; x2], [d1], d_rec});
part.add_mode({fx{3}, [x1; x2], [d1], d_rec});
part.add_mode({fx{4}, [x1; x2], [d1], d_rec});

clf
plot_vf(part)

