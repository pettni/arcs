u_coolant_valve = 1;
u_radiator_grill = 1;
Q_engine_combustion_heat = 2E4; %[1E4, 2E4]; %% [W]
Ta = 270; %[260 300]; %% [K]
W_coolant_pump = 0.02; %% [0.02, 0.05];  %% [kg/s]
vehicle_speed_mps = 20; %[10, 20];  %% [m/s]

cp_air = 1005; %% J/kg/K   (specific heat of air [1005])
cp_coolant = 3400; %% J/kg/K   (specific heat of coolant [3400])

C_engine_block = 750;  %% J/K (Assume 60% AL and 40% FE)
C_radiator = 200; %% J/K (Assume all aluminium)

U_block2amb = 100;    %% W/K (Tuned to match responses)%%
U_radiator2amb = 100; %% W/K (Tuned to match responses)   

radiator_frontal_area = 0.2; %% m^2
density_air = 1; %%kg/m2

W_airflow_front = vehicle_speed_mps * radiator_frontal_area * density_air; %%kg/s
W_airflow_radiator = W_airflow_front * u_radiator_grill;
W_coolant_radiator = W_coolant_pump*u_coolant_valve;

Ao11 = [[-U_block2amb - cp_coolant*W_coolant_radiator,  cp_coolant*W_coolant_radiator]/C_engine_block;
    [cp_coolant * W_coolant_radiator, -cp_air*W_airflow_radiator - cp_coolant*W_coolant_radiator- U_radiator2amb]/C_radiator];
Ko11 = [(Q_engine_combustion_heat+U_block2amb*Ta)/C_engine_block;(cp_air*W_airflow_radiator+U_radiator2amb)*Ta/C_radiator];


W_airflow_radiator = W_airflow_front * 0.25;
W_coolant_radiator = W_coolant_pump*1;

Ao12 = [[-U_block2amb - cp_coolant*W_coolant_radiator,  cp_coolant*W_coolant_radiator]/C_engine_block;
    [cp_coolant * W_coolant_radiator, -cp_air*W_airflow_radiator - cp_coolant*W_coolant_radiator- U_radiator2amb]/C_radiator];
Ko12 = [(Q_engine_combustion_heat+U_block2amb*Ta)/C_engine_block;(cp_air*W_airflow_radiator+U_radiator2amb)*Ta/C_radiator];


W_airflow_radiator = W_airflow_front *1;
W_coolant_radiator = W_coolant_pump*0.25;

Ao21 = [[-U_block2amb - cp_coolant*W_coolant_radiator,  cp_coolant*W_coolant_radiator]/C_engine_block;
    [cp_coolant * W_coolant_radiator, -cp_air*W_airflow_radiator - cp_coolant*W_coolant_radiator- U_radiator2amb]/C_radiator];
Ko21 = [(Q_engine_combustion_heat+U_block2amb*Ta)/C_engine_block;(cp_air*W_airflow_radiator+U_radiator2amb)*Ta/C_radiator];


W_airflow_radiator = W_airflow_front * 0.25;
W_coolant_radiator = W_coolant_pump*0.25;

Ao22 = [[-U_block2amb - cp_coolant*W_coolant_radiator,  cp_coolant*W_coolant_radiator]/C_engine_block;
    [cp_coolant * W_coolant_radiator, -cp_air*W_airflow_radiator - cp_coolant*W_coolant_radiator- U_radiator2amb]/C_radiator];
Ko22 = [(Q_engine_combustion_heat+U_block2amb*Ta)/C_engine_block;(cp_air*W_airflow_radiator+U_radiator2amb)*Ta/C_radiator];

domain = Rec([250 250; 
              500 500]);
goal_set = Rec([380 250; 
                390 500], 1);
unsafe_set = Rec([400 250; 
                  500 500], 2);
% Model
fx1.A = Ao11; fx1.K = Ko11;
fx2.A = Ao12; fx2.K = Ko12;
fx3.A = Ao21; fx3.K = Ko21;
fx4.A = Ao22; fx4.K = Ko22;

% Action set
act_set={fx1,fx2,fx3,fx4};

% Build initial partition
part = Partition(domain);
part.add_area(goal_set);
part.add_area(unsafe_set)
part.check();   % sanity check