%% Main function to generate tests
function tests = exampleTest
  tests = functiontests(localfunctions);
end

function test_lone_build1(testCase)
  pos = [0,0];
  slabs = [0];
  room_types = [1];
  
  R_o = 2.356*18/36 / 24;
  C = 5*1012*1.2041*(8*6*3);
  T_o = 30;
  
  answer_A = -1/(C*R_o);
  answer_K = 1/(C*R_o)*T_o + 6*(8*6)/C;
  answer_E = 1/C*(8*6);
  
  build = Building.create_building(pos, slabs, room_types);
  [A, E, K] = build.get_dyn();
  
  verifyEqual(testCase, A{1}, answer_A, 'RelTol', 0.001);
  verifyEqual(testCase, E{1}, answer_E, 'RelTol', 0.001);
  verifyEqual(testCase, K{1}, answer_K, 'RelTol', 0.001);
  clear all
end

function test_lone_build2(testCase)
  pos = [0,0];
  slabs = [0];
  room_types = [2];
  
  R_o = 2.356*18/36 / 18;
  C = 5*1012*1.2041*(6*6*3);
  T_o = 30;
  
  answer_A = -1/(C*R_o);
  answer_K = 1/(C*R_o)*T_o + 8*(6*6)/C;
  answer_E = 1/C*(6*6);
  
  build = Building.create_building(pos, slabs, room_types);
  [A, E, K] = build.get_dyn();
  
  verifyEqual(testCase, A{1}, answer_A, 'RelTol', 0.001);
  verifyEqual(testCase, E{1}, answer_E, 'RelTol', 0.001);
  verifyEqual(testCase, K{1}, answer_K, 'RelTol', 0.001);
  clear all
end

function test_lone_build_slabs(testCase)    
    floor_A = [8*6, 6*6];
    A_fac = [24, 18];
    V = floor_A*3;
    heats = [6, 8];
    
    for id = 1:2
        pos = [0,0];
        slabs = [1];
        room_types = [id];
        
        C_c = 0.25*2300*750*floor_A(id);
        R_c = 0.125/floor_A(id);
        R_w = 0.102/floor_A(id);
        R_o = 2.356 * 18/36 / A_fac(id);
        C_1 = 5*1012*1.2041*V(id);
        T_w = 18;
        T_o = 30;
        
        answer_A = @(R) [-(1/(C_c*R_c) + 1/(C_c*R)), 1/(C_c*R_c);
                    1/(C_1*R_c), -(1/(C_1*R_c) + 1/(C_1*R_o))];
        answer_K = @(R) [1/(C_c*R)*T_w; heats(id)*floor_A(id)/C_1 + 1/(C_1*R_o)*T_o];
        
        answer_A1 = answer_A(R_w);
        answer_A2 = answer_A(Inf);
        answer_K1 = answer_K(R_w);
        answer_K2 = answer_K(Inf);
        
        answer_E = [0; floor_A(id)/C_1];
        
        clear build;
        build = Building.create_building(pos, slabs, room_types);
        [A,E,K] = build.get_dyn();
        
        verifyEqual(testCase, A, {answer_A1, answer_A2}, 'RelTol', 0.00001);
        verifyEqual(testCase, E, {answer_E, answer_E}, 'RelTol', 0.00001);
        verifyEqual(testCase, K, {answer_K1, answer_K2}, 'RelTol', 0.00001);
    end
    clear all
end

function test_building_same_slab(testCase)
    pos = [0,0;
           1,0];
    slabs = [1, 1];
    room_types = [1,2];
    
    answer_A1 = [-0.0413,    0.0106,    0.0080;
                0.4377,   -0.4869,    0.0260;
                0.4377,    0.0346,   -0.4955]*1e-3;
    answer_E1 = [0,         0;
                0.5471,         0;
                0,    0.5471]*1e-4;
    answer_K1 = [0.0004;
                0.0010;
                0.0011];
    answer_A2 = [-0.0186,    0.0106,    0.0080;
    0.4377,   -0.4869,    0.0260;
    0.4377,    0.0346,   -0.4955]*1e-3;

    answer_E2 = [0,         0;
    0.5471,         0;
         0,    0.5471]*1e-4;
    
    answer_K2 = [0;
    0.0010;
    0.0011];
            
    build = Building.create_building(pos, slabs, room_types);
    [A, E, K] = build.get_dyn();
    
    verifyEqual(testCase, A, {answer_A1, answer_A2}, 'RelTol', 0.05);
    verifyEqual(testCase, E, {answer_E1, answer_E2}, 'RelTol', 0.05);
    verifyEqual(testCase, K, {answer_K1, answer_K2}, 'RelTol', 0.05);
end

function test_rooms_only_one_slab(testCase)
    pos = [0,0;
           1,0];
    slabs = [1,0];
    room_types = [1,2];
    
    A_f1 = 8*6;
    A_f2 = 6*6;
    A_fac1 = 24;
    A_fac2 = 18;
    
    C_c = 0.25*750*2300*A_f1;
    R_c = 0.125/A_f1;
    R_w = 0.102/A_f1;
    C_1 = 5*1012*1.2041*A_f1*3;
    C_2 = 5*1012*1.2041*A_f2*3;
    R_12 = 0.79/18;
    R_o1 = 2.356*18/36 / A_fac1;
    R_o2 = 2.356*18/36 / A_fac2;
    q_1 = 6*A_f1;
    q_2 = 8*A_f2;
    T_w = 18;
    T_o = 30;
    
    answer_A = @(R) [-(1/(C_c*R_c) + 1/(C_c*R)), 1/(R_c*C_c), 0;
                 1/(C_1*R_c), -(1/(C_1*R_12) + 1/(C_1*R_o1) + 1/(C_1*R_c)), 1/(C_1*R_12);
                 0, 1/(C_2*R_12), -(1/(C_2*R_12) + 1/(C_2*R_o2))]; 
    
    answer_K = @(R) [T_w/(R*C_c);
                 q_1/C_1 + T_o/(C_1*R_o1);
                 q_2/C_2 + T_o/(C_2*R_o2)];
             
    answer_A1 = answer_A(R_w);
             
    answer_A2 = answer_A(Inf);
    
    answer_K1 = answer_K(R_w);
    
    answer_K2 = answer_K(Inf);
    
    answer_E = [0, 0;
                A_f1/C_1, 0;
                0, A_f2/C_2];
             
    build = Building.create_building(pos, slabs, room_types);
    [A, E, K] = build.get_dyn();
    
    verifyEqual(testCase, A, {answer_A1, answer_A2}, 'RelTol', 0.00001);
    verifyEqual(testCase, K, {answer_K1, answer_K2}, 'RelTol', 0.00001);
    verifyEqual(testCase, E, {answer_E, answer_E}, 'RelTol', 0.00001);
   
end

function test_two_room_slabs(testCase)
    pos = [0,0;
           1,0];
    slabs = [1,2];
    room_types = [1,2];
    
    A_f1 = 8*6;
    A_f2 = 6*6;
    A_fac1 = 24;
    A_fac2 = 18;
    
    C_c1 = 0.25*750*2300*A_f1;
    C_c2 = 0.25*750*2300*A_f2;
    R_c1 = 0.125/A_f1;
    R_c2 = 0.125/A_f2;
    R_w1 = 0.102/A_f1;
    R_w2 = 0.102/A_f2;
    C_1 = 5*1012*1.2041*A_f1*3;
    C_2 = 5*1012*1.2041*A_f2*3;
    R_12 = 0.79/18;
    R_o1 = 2.356*18/36 / A_fac1;
    R_o2 = 2.356*18/36 / A_fac2;
    q_1 = 6*A_f1;
    q_2 = 8*A_f2;
    T_w = 18;
    T_o = 30;
    
    answer_A = @(R1, R2) [-1/C_c1*(1/R_c1 + 1/R1), 0, 1/(C_c1*R_c1), 0;
                          0, -1/C_c2*(1/R_c2 + 1/R2), 0, 1/(C_c2*R_c2);
                          1/(C_1*R_c1), 0, -1/C_1*(1/R_c1 + 1/R_12 + 1/R_o1), 1/(C_1*R_12);
                          0, 1/(C_2*R_c2), 1/(C_2*R_12), -1/C_2*(1/R_c2 + 1/R_12 + 1/R_o2)];
    
    answer_K = @(R1, R2) [T_w/(C_c1*R1); T_w/(C_c2*R2);
                          q_1/C_1 + T_o/(C_1*R_o1); q_2/C_2 + T_o/(C_2*R_o2)];
    
    answer_A1 = answer_A(R_w1, R_w2);
    answer_A2 = answer_A(R_w1, Inf);
    answer_A3 = answer_A(Inf, R_w2);
    answer_A4 = answer_A(Inf, Inf);
    
    answer_K1 = answer_K(R_w1, R_w2);
    answer_K2 = answer_K(R_w1, Inf);
    answer_K3 = answer_K(Inf, R_w2);
    answer_K4 = answer_K(Inf, Inf);
    
    answer_E = [0, 0;
                0, 0;
                A_f1/C_1, 0;
                0, A_f2/C_2];
    
    build = Building.create_building(pos, slabs, room_types);
    [A, E, K] = build.get_dyn();
    
    verifyEqual(testCase, A, {answer_A1,answer_A2,answer_A3,answer_A4}, 'RelTol', 0.00001);
    verifyEqual(testCase, K, {answer_K1,answer_K2,answer_K3,answer_K4}, 'RelTol', 0.00001);
    verifyEqual(testCase, E, {answer_E, answer_E, answer_E, answer_E}, 'RelTol', 0.00001);
end

function test_multiroom_slab(testCase)
    pos = [0,0;
           1,0;
           1,1];
    slabs = [1,1,1];
    room_types = [1,2,1];
    
    A_f1 = 8*6;
    A_f2 = 6*6;
    A_fac1 = 24;
    A_fac2 = 18;
    
    C_c = 0.25*750*2300*(2*A_f1 + A_f2);
    R_c1 = 0.125/A_f1;
    R_c2 = 0.125/A_f2;
    R_c3 = R_c1;
    R_w = 0.102/(2*A_f1+A_f2);
    C_1 = 5*1012*1.2041*A_f1*3;
    C_2 = 5*1012*1.2041*A_f2*3;
    C_3 = C_1;
    R_12 = 0.79/18;
    R_23 = R_12;
    R_o1 = 2.356*18/36 / A_fac1;
    R_o2 = 2.356*18/36 / A_fac2;
    R_o3 = R_o1;
    q_1 = 6*A_f1;
    q_2 = 8*A_f2;
    q_3 = q_1;
    T_w = 18;
    T_o = 30;
    
    answer_A = @(R) [-1/C_c*(1/R_c1 + 1/R_c2 + 1/R_c3 + 1/R), 1/(C_c*R_c1), 1/(C_c*R_c2), 1/(C_c*R_c3);
        1/(C_1*R_c1), -1/C_1*(1/R_c1 + 1/R_12 + 1/R_o1), 1/(C_1*R_12), 0;
        1/(C_2*R_c2), 1/(C_2*R_12), -1/C_2*(1/R_c2 + 1/R_12 + 1/R_23 + 1/R_o2), 1/(C_2*R_23);
        1/(C_3*R_c1), 0, 1/(C_3*R_23), -1/C_3*(1/R_c1 + 1/R_23 + 1/R_o3)];
    
    answer_K = @(R) [T_w/(C_c*R); T_o/(C_1*R_o1) + q_1/C_1; T_o/(C_2*R_o2) + q_2/C_2; T_o/(C_3*R_o3) + q_3/C_3];
    
    
    answer_A1 = answer_A(R_w);
    answer_A2 = answer_A(Inf);
    
    answer_K1 = answer_K(R_w);
    answer_K2 = answer_K(Inf);
    
    answer_E = [0, 0, 0;
                A_f1/C_1, 0, 0;
                0, A_f2/C_2, 0;
                0, 0, A_f1/C_3];
    
    build = Building.create_building(pos, slabs, room_types);
    [A, E, K] = build.get_dyn();
    
    verifyEqual(testCase, A, {answer_A1,answer_A2}, 'RelTol', 0.00001);
    verifyEqual(testCase, K, {answer_K1,answer_K2}, 'RelTol', 0.00001);
    verifyEqual(testCase, E, {answer_E, answer_E}, 'RelTol', 0.00001);
end