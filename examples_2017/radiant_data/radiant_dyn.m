%% get_dyn: get TABS dynamics
function [DYN_A1, DYN_K1, DYN_E1, DYN_A2, DYN_K2, DYN_E2] = radiant_dyn()

	% Room 1 parameters
	d_x1 = 8;			  % [m] room width
	d_y1 = 6;			  % [m] room depth
	d_z1 = 3;				% [m] room height
	A_fac1 = 24;		% [m^2] facade area
	R_cov1 = 0.125; % [m^2K/W] floor thermal resistance
	n_a1 = 0.1/3600;% [1/s] floor infiltration coefficient

	% Room 2 parameters
	d_x2 = 6;			  % [m] room width
	d_y2 = 6;			  % [m] room depth
	d_z2 = 3;				% [m] room height
	A_fac2 = 18; 		% [m^2] facade area
	R_cov2 = 0.125;  % [m^2K/W] floor thermal resistance
	n_a2 = 0.1/3600;% [1/s] floor infiltration coefficient

	A_12 = 18; 			% [m^2] area of wall connecting rooms
	U = 0.65;				% [W/m^2K] facade U value

	% Piping/slab system
	d_s = 0.25; 		% [m] thickness of slab

	% Natural constants
	rho_c = 2300;		% [kg/m3] density, concrete
	c_c = 750;			% [J/kgK] specific heat, concrete
	rho_a = 1.2041;	% [kg/m^3] density, air
	c_a = 1012; 		% [J/kgK] specific heat, air		
	lambda_s = 1;   % [W/mK] thermal conductivity lightweight concrete [0.1-2] for different densities

	A_fl1 = d_x1 * d_y1;
	V_1 = A_fl1 * d_z1;
	A_fl2 = d_x2 * d_y2;
	V_2 = A_fl2 * d_z2;
	A_slab = A_fl1 + A_fl2;

	% These formulas are from M. Gwerder et al. / Applied Energy 85 (2008) 565–581 
	% Room/slab R values [m^2K/W]
	R_tilde = 0.125;
	% Room/outside R value [m^2K/W]  
	% (adjusted: now per unit facade area instead of unit floor area)
	R_lf = 2.356*18/36;
	% Water/slab R value
	R_t = 0.102;  % [m^2K/W]

	% Room/room R value
	% from "Building models for model predictive control of office buildings with concrete core activation"
	r_12 = 0.79;  % [m^2K/W]

	% Thermal heat insulation [m^K/W]
	r_1o = R_lf;
	r_2o = R_lf;
	r_1c = R_tilde;
	r_2c = R_tilde;
	r_cw = R_t;

	% Thermal resistance [K/W]
	R_1o = r_1o/A_fac1;
	R_2o = r_2o/A_fac2;
	R_12 = r_12/A_12;
	R_1c = r_1c/A_fl1;
	R_2c = r_2c/A_fl2;
	R_cw = r_cw/A_slab;

	K1 = 1/R_1o;
	K2 = 1/R_2o;
	K12 = 1/R_12;
	K21 = 1/R_12;
	Kw = 1/R_cw;
	Kr1 = 1/R_1c;
	Kr2 = 1/R_2c;

	% Thermal capacitances [J/K]
	C1 = 5 * rho_a * c_a * V_1;
	C2 = 5 * rho_a * c_a * V_2;
	Cr = rho_c * c_c * d_s * (A_slab);

	% Temperatures [C]
	Ta = 30;
	Tw = 18;

	% Added heat [J/m^2]
	q1 = 6;
	q2 = 8;

	disp(['$r_{1o}=r_{2o}=', num2str(r_1o), 'm^2K/W$, ', ...
			  '$r_{1c}=r_{2c}=', num2str(r_1c), 'm^2K/W$, ', ...
			  '$r_{cw}=', num2str(r_cw), 'm^2K/W$'])
	disp(['$R_{1c} = \frac{r_{1c}}{A_{1c}} = ', num2str(R_1c), 'K/W$'])
	disp(['$R_{1o} = \frac{r_{1o}}{A_{1o}} = ', num2str(R_1o), 'K/W$'])
	disp(['$R_{2c} = \frac{r_{2c}}{A_{wc}} = ', num2str(R_2c), 'K/W$'])
	disp(['$R_{2o} = \frac{r_{2o}}{A_{2o}} = ', num2str(R_2o), 'K/W$'])
	disp(['$R_{12} = \frac{r_{12}}{A_{12}} = ', num2str(R_12), 'K/W$'])
	disp(['$R_{cw} = \frac{r_{cw}}{A_{cw}} = ', num2str(R_cw), 'K/W$'])
	disp(['$C_{1} = 5 c_a \rho_a V_1 = ', num2str(C1), 'J/K$'])
	disp(['$C_{2} = 5 c_a \rho_a V_2 = ', num2str(C2), 'J/K$'])
	disp(['$C_{c} = c_c \rho_c V_c = ', num2str(Cr), 'J/K$'])
	disp(['$q_1 = \left(', num2str(q1) , '\frac{W}{m^2} \right) \times A_{1c}$'])
	disp(['$q_2 = \left(', num2str(q2) , '\frac{W}{m^2} \right) \times A_{2c}$'])

	DYN_A1 = diag([1/Cr, 1/C1, 1/C2]) * ...
					 [-(Kr1+Kr2+Kw)	Kr1 					Kr2;	
				  	Kr1 			-(Kr1+K1+K12)	K12;
						Kr2				K21						-(Kr2+K2+K21)];

	DYN_K1 = diag([1/Cr, 1/C1, 1/C2]) * ...
				   [Kw*Tw;
					  K1*Ta+A_fl1*q1;
					  K2*Ta+A_fl2*q2];

	DYN_E1 = diag([1/Cr, 1/C1, 1/C2]) * ...
					 [0 0;
					 	A_fl1 0;
					 	0 A_fl2];

	DYN_A2 = diag([1/Cr, 1/C1, 1/C2]) * ...
						 [-(Kr1+Kr2)  	Kr1       		Kr2;  
	    		  	Kr1     	   -(Kr1+K1+K12)  K12;
	    				Kr2 		    	K21      			-(Kr2+K2+K21)];

	DYN_K2 = diag([1/Cr, 1/C1, 1/C2]) * ...
					 [0;
		    		K1*Ta+A_fl1*q1;
		    		K2*Ta+A_fl2*q2];

	DYN_E2 = diag([1/Cr, 1/C1, 1/C2]) * ...
					 [0 0;
					 	A_fl1 0;
					 	0 A_fl2];

end