domain = Rec([20 28; 20 28; 20 28]);
goal_set = Rec([20, 28; 22 26; 22 26], 1);
unsafe_set = Rec([27.7 28; 27 28; 27.2 28], 2);

load('radiant_data/a1.mat')
load('radiant_data/a2.mat')
load('radiant_data/b1.mat')
load('radiant_data/b2.mat')

b1(3) = b1(3)/10; 
b2(3) = b2(3)/10;

fx1.A = a1;
fx1.K = b1;

fx2.A = a2;
fx2.K = b2;
act_set={fx1,fx2};

% Build initial partition
part = Partition(domain);
part.add_area(goal_set);
part.add_area(unsafe_set)
part.check();   % sanity check