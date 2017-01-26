%% get_dyn: get TABS dynamics
function [DYN_A1, DYN_K1, DYN_A2, DYN_K2] = radiant_dyn()
	% Room 1 parameters
	d_x1 = 6;			  % [m] room width
	d_y1 = 6;			  % [m] room depth
	d_z1 = 3;				% [m] room height
	A_fac1 = 18;		% [m^2] facade area
	R_cov1 = 0.125; % [m^2K/W] floor thermal resistance
	n_a1 = 0.1/3600;% [1/s] floor infiltration coefficient

	% Room 2 parameters
	d_x2 = 6;			  % [m] room width
	d_y2 = 6;			  % [m] room depth
	d_z2 = 3;				% [m] room height
	A_fac2 = 18; 		% [m^2] facade area
	R_cov2 = 0.125; % [m^2K/W] floor thermal resistance
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
	c_w = 4186;			% [J/kgK] specific heat water
	lambda_s = 1;   % [W/mK] thermal conductivity lightweight concrete [0.1-2] for different densities

	A_fl1 = d_x1 * d_y1;
	V_1 = A_fl1 * d_z1;
	A_fl2 = d_x2 * d_y2;
	V_2 = A_fl2 * d_z2;
	A_slab = A_fl1 + A_fl2;

	% These formulas are from M. Gwerder et al. / Applied Energy 85 (2008) 565–581 
	% Room/slab R values [m^2K/W]
	R_tilde1 = ((d_s/2)/lambda_s + R_cov1)/2;
	R_tilde2 = ((d_s/2)/lambda_s + R_cov2)/2;

	% Room/outside R value [m^2K/W]  (per unit facade area)
	R_lf1 = 1/(U + n_a1 * V_1 * rho_a * c_a / A_fac1);
	R_lf2 = 1/(U + n_a2 * V_2 * rho_a * c_a / A_fac2);

	% Water/slab R value (formula is complicated)
	R_t = 0.102;  % [m^2K/W]

	% Room/room R value
	% from "Building models for model predictive control of office buildings with concrete core activation"
	R_12 = 0.79;  % [m^2K/W]

	% Temperatures [C]
	Ta = 30;
	Tw = 18;

	% Thermal conductances [W/K]
	K1 = A_fac1/R_lf1;
	K2 = A_fac2/R_lf2;
	K12 = A_12/R_12;
	K21 = A_12/R_12;
	Kw = A_slab/R_t;
	Kr1 = A_fl1/R_tilde1;
	Kr2 = A_fl2/R_tilde2;

	% Thermal capacitances [J/K]
	C1 = 5 * rho_a * c_a * V_1;
	C2 = 5 * rho_a * c_a * V_2;
	Cr = rho_c * c_c * d_s * (A_fl1 + A_fl2);

	% Added heat [J]
	q1 = 6 * A_fl1;
	q2 = 8 * A_fl2;

	DYN_A1 = diag([1/Cr, 1/C1, 1/C2]) * ...
					 [-(Kr1+Kr2+Kw)	Kr1 					Kr2;	
				  	Kr1 			-(Kr1+K1+K12)	K12;
						Kr2				K21						-(Kr2+K2+K21)];

	DYN_K1 = diag([1/Cr, 1/C1, 1/C2]) * ...
				   [Kw*Tw;
					  K1*Ta+q1;
					  K2*Ta+q2];

	DYN_A2 = diag([1/Cr, 1/C1, 1/C2]) * ...
						 [-(Kr1+Kr2)  	Kr1       		Kr2;  
	    		  	Kr1     	   -(Kr1+K1+K12)  K12;
	    				Kr2 		    	K21      			-(Kr2+K2+K21)];

	DYN_K2 = diag([1/Cr, 1/C1, 1/C2]) * ...
					 [0;
		    		K1*Ta+q1;
		    		K2*Ta+q2];
end