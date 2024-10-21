function [A, b, center_1, center_2, center_3] = Matrix_definition(nacapane_1, nacapane_2, nacapane_3, N1, N2, N3, u_inf)


center_1  = zeros(N1,2);                                                     % MEMORY ALLOCATION
center_2  = zeros(N2,2);                                                     % MEMORY ALLOCATION
center_3  = zeros(N3,2);                                                     % MEMORY ALLOCATION
A_s_11     = zeros(N1);                                                      % MEMORY ALLOCATION
A_s_12     = zeros(N1, N2);                                                  % MEMORY ALLOCATION
A_s_13     = zeros(N1, N3);                                                  % MEMORY ALLOCATION
A_s_21     = zeros(N2, N1);                                                  % MEMORY ALLOCATION
A_s_22     = zeros(N2);                                                      % MEMORY ALLOCATION
A_s_23     = zeros(N2, N3);                                                  % MEMORY ALLOCATION
A_s_31     = zeros(N3, N1);                                                  % MEMORY ALLOCATION
A_s_32     = zeros(N3, N2);                                                  % MEMORY ALLOCATION
A_s_33     = zeros(N3);                                                      % MEMORY ALLOCATION

b_s_1     = zeros(N1, 1);                                                    % MEMORY ALLOCATION
b_s_2     = zeros(N2, 1);                                                    % MEMORY ALLOCATION
b_s_3     = zeros(N3, 1);                                                    % MEMORY ALLOCATION

a_v_11     = zeros(N1, 1);                                                   % MEMORY ALLOCATION
a_v_12     = zeros(N1, 1);                                                   % MEMORY ALLOCATION
a_v_13     = zeros(N1, 1);                                                   % MEMORY ALLOCATION
a_v_21     = zeros(N2, 1);                                                   % MEMORY ALLOCATION
a_v_22     = zeros(N2, 1);                                                   % MEMORY ALLOCATION
a_v_23     = zeros(N2, 1);                                                   % MEMORY ALLOCATION
a_v_31     = zeros(N3, 1);                                                   % MEMORY ALLOCATION
a_v_32     = zeros(N3, 1);                                                   % MEMORY ALLOCATION
a_v_33     = zeros(N3, 1);                                                   % MEMORY ALLOCATION


% POINTS FOR KUTTA CONDITION - FOIL 1
point_1_1 = [nacapane_1.x_pane(1); nacapane_1.y_pane(1)];                   % DEFINITION OF COORDINATES OF FIRST POINT OF FIRST PANEL FOIL 1
point_2_1 = [nacapane_1.x_pane(2); nacapane_1.y_pane(2)];                   % DEFINITION OF COORDINATES OF LAST POINT OF FIRST PANEL FOIL 1
point_N_1 = [nacapane_1.x_pane(N1); nacapane_1.y_pane(N1)];                 % DEFINITION OF COORDINATES OF FIRST POINT OF LAST PANEL FOIL 1
point_NN_1 = [nacapane_1.x_pane(N1+1); nacapane_1.y_pane(N1+1)];            % DEFINITION OF COORDINATES OF LAST POINT OF LAST PANEL FOIL 1

[~, tau1_1] = versors(point_1_1, point_2_1);                                % DEFINITION OF TANGENT VERSOR OF FIRST PANEL FOIL 1
[~, tauN_1] = versors(point_N_1, point_NN_1);                               % DEFINITION OF TANGENT VERSOR OF LAST PANEL FOIL 1

% POINTS FOR KUTTA CONDITION - FOIL 2
point_1_2 = [nacapane_2.x_pane(1); nacapane_2.y_pane(1)];                   % DEFINITION OF COORDINATES OF FIRST POINT OF FIRST PANEL FOIL 2
point_2_2 = [nacapane_2.x_pane(2); nacapane_2.y_pane(2)];                   % DEFINITION OF COORDINATES OF LAST POINT OF FIRST PANEL FOIL 2
point_N_2 = [nacapane_2.x_pane(N2); nacapane_2.y_pane(N2)];                 % DEFINITION OF COORDINATES OF FIRST POINT OF LAST PANEL FOIL 2
point_NN_2 = [nacapane_2.x_pane(N2+1); nacapane_2.y_pane(N2+1)];            % DEFINITION OF COORDINATES OF LAST POINT OF LAST PANEL FOIL 2

[~, tau1_2] = versors(point_1_2, point_2_2);                                % DEFINITION OF TANGENT VERSOR OF FIRST PANEL FOIL 2
[~, tauN_2] = versors(point_N_2, point_NN_2);                               % DEFINITION OF TANGENT VERSOR OF LAST PANEL FOIL 2

% POINTS FOR KUTTA CONDITION - FOIL 3
point_1_3 = [nacapane_3.x_pane(1); nacapane_3.y_pane(1)];                   % DEFINITION OF COORDINATES OF FIRST POINT OF FIRST PANEL FOIL 2
point_2_3 = [nacapane_3.x_pane(2); nacapane_3.y_pane(2)];                   % DEFINITION OF COORDINATES OF LAST POINT OF FIRST PANEL FOIL 2
point_N_3 = [nacapane_3.x_pane(N3); nacapane_3.y_pane(N3)];                 % DEFINITION OF COORDINATES OF FIRST POINT OF LAST PANEL FOIL 2
point_NN_3 = [nacapane_3.x_pane(N3+1); nacapane_3.y_pane(N3+1)];            % DEFINITION OF COORDINATES OF LAST POINT OF LAST PANEL FOIL 2

[~, tau1_3] = versors(point_1_3, point_2_3);                                % DEFINITION OF TANGENT VERSOR OF FIRST PANEL FOIL 2
[~, tauN_3] = versors(point_N_3, point_NN_3);                               % DEFINITION OF TANGENT VERSOR OF LAST PANEL FOIL 2

% 11

c_s_1 = zeros(1,N1);
c_s_N = zeros(1,N1);
c_v_1 = zeros(1,N1);
c_v_N = zeros(1,N1);

for i = 1 : N1
    point_i = [nacapane_1.x_pane(i); nacapane_1.y_pane(i)];                 % DEFINITION OF COORDINATES OF FIRST POINT OF PANEL i
    point_ii = [nacapane_1.x_pane(i+1); nacapane_1.y_pane(i+1)];            % DEFINITION OF COORDINATES OF LAST POINT OF PANEL i
    [n, ~] = versors (point_i, point_ii);                                   % DEFINITION OF NORMAL VERSOR OF PANEL i
    center_x = (point_i(1) + point_ii(1))*0.5;                              % DEFINITION OF X-COORDINATE OF CENTER POINT OF PANEL i
    center_y = (point_i(2) + point_ii(2))*0.5;                              % DEFINITION OF X-COORDINATE OF CENTER POINT OF PANEL i
    Centro = [center_x; center_y];                                          % DEFINITION OF COORDINATES OF CENTER POINT OF PANEL i
    center_1(i, :) = Centro;                                                % STORING THE COORDINATE OF CENTER POINTS OF ALL PANELS
    
    for j = 1 : N1
        Estremo_1 = [nacapane_1.x_pane(j); nacapane_1.y_pane(j)];           % DEFINITION OF COORDINATES OF FIRST POINT OF PANEL j               
        Estremo_2 = [nacapane_1.x_pane(j+1); nacapane_1.y_pane(j+1)];       % DEFINITION OF COORDINATES OF LAST POINT OF PANEL j
        beta = atan2((Estremo_2(2) - Estremo_1(2)) , ...
                     (Estremo_2(1) - Estremo_1(1)));                        % LOCAL SDR ANGLE [rad]
        Q = [ cos(beta), sin(beta);
             -sin(beta), cos(beta)                                          % G2L MATRIX
            ];                                                              

        Estremo_1_m = [nacapane_1.x_pane_m(j); nacapane_1.y_pane_m(j)];         % DEFINITION OF COORDINATES OF FIRST POINT OF PANEL j               
        Estremo_2_m = [nacapane_1.x_pane_m(j+1); nacapane_1.y_pane_m(j+1)];     % DEFINITION OF COORDINATES OF LAST POINT OF PANEL j
        beta_m = atan2((Estremo_2_m(2) - Estremo_1_m(2)) , ...
                     (Estremo_2_m(1) - Estremo_1_m(1)));                    % LOCAL SDR ANGLE [rad]
        Q_m = [ cos(beta_m), sin(beta_m);
             -sin(beta_m), cos(beta_m)                                      % G2L MATRIX
            ]; 

        u_sorg = ViSorgente(Centro, Estremo_1, Estremo_2, Q', Q);           % COMPUTATION OF VELOCITY ON CENTER OF PANEL i DUE TO PANEL SOURCE j
        u_sorg_m = ViSorgente(Centro, Estremo_1_m, Estremo_2_m, Q_m', Q_m); % COMPUTATION OF VELOCITY ON CENTER OF PANEL i DUE TO MIRROR PANEL SOURCE j        
        u_vort = ViVortice(Centro, Estremo_1, Estremo_2, Q', Q);            % COMPUTATION OF VELOCITY ON CENTER OF PANEL i DUE TO PANEL VORTEX j
        u_vort_m = ViVortice(Centro, Estremo_1_m, Estremo_2_m, Q_m', Q_m);  % COMPUTATION OF VELOCITY ON CENTER OF PANEL i DUE TO MIRROR PANEL VORTEX j

        A_s_11(i,j) = (u_sorg' + u_sorg_m') * n;                                          % DEFINITION OF ELEMENT i,j OF A_S SUBMATRIX
        a_v_11(i) = a_v_11(i) + (u_vort' - u_vort_m') * n;                                     % DEFINITION OF ELEMENT i OF a_v 
        
        if i == 1                                                           % IF-CYCLE TO DEFINE c_s AND c_v COEFFICIENT
          c_s_1(j) = (u_sorg' + u_sorg_m') * tau1_1;
          c_v_1(j) = (u_vort' - u_vort_m') * tau1_1;
        end
        if i == N1
          c_s_N(j) = (u_sorg' + u_sorg_m') * tauN_1;
          c_v_N(j) = (u_vort' - u_vort_m') * tauN_1;
        end
    end

    b_s_1(i) = -u_inf' * n;                                                 % DEFINITION OF b_s VECTOR
    b_v_1 = - u_inf' * (tau1_1 + tauN_1);                                   % DEFINITION OF b_v VARIABLE

end
c_s_11 = c_s_1 + c_s_N;                                                     % FINAL CHARACTERIZATION OF c_s
c_v_11 = sum(c_v_1) + sum(c_v_N);                                           % FINAL CHARACTERIZATION OF c_s

% 12

c_s_1 = zeros(1,N2);
c_s_N = zeros(1,N2);
c_v_1 = zeros(1,N2);
c_v_N = zeros(1,N2);

for i = 1 : N1
    point_i = [nacapane_1.x_pane(i); nacapane_1.y_pane(i)];                 % DEFINITION OF COORDINATES OF FIRST POINT OF PANEL i
    point_ii = [nacapane_1.x_pane(i+1); nacapane_1.y_pane(i+1)];            % DEFINITION OF COORDINATES OF LAST POINT OF PANEL i
    [n, ~] = versors (point_i, point_ii);                                   % DEFINITION OF NORMAL VERSOR OF PANEL i
    center_x = (point_i(1) + point_ii(1))*0.5;                              % DEFINITION OF X-COORDINATE OF CENTER POINT OF PANEL i
    center_y = (point_i(2) + point_ii(2))*0.5;                              % DEFINITION OF X-COORDINATE OF CENTER POINT OF PANEL i
    Centro = [center_x; center_y];                                          % DEFINITION OF COORDINATES OF CENTER POINT OF PANEL i
    
    for j = 1 : N2
        Estremo_1 = [nacapane_2.x_pane(j); nacapane_2.y_pane(j)];           % DEFINITION OF COORDINATES OF FIRST POINT OF PANEL j               
        Estremo_2 = [nacapane_2.x_pane(j+1); nacapane_2.y_pane(j+1)];       % DEFINITION OF COORDINATES OF LAST POINT OF PANEL j
        beta = atan2((Estremo_2(2) - Estremo_1(2)) , ...
                     (Estremo_2(1) - Estremo_1(1)));                        % LOCAL SDR ANGLE [rad]
        Q = [ cos(beta), sin(beta);
             -sin(beta), cos(beta)                                          % G2L MATRIX
            ];                

        Estremo_1_m = [nacapane_2.x_pane_m(j); nacapane_2.y_pane_m(j)];         % DEFINITION OF COORDINATES OF FIRST POINT OF PANEL j               
        Estremo_2_m = [nacapane_2.x_pane_m(j+1); nacapane_2.y_pane_m(j+1)];     % DEFINITION OF COORDINATES OF LAST POINT OF PANEL j
        beta_m = atan2((Estremo_2_m(2) - Estremo_1_m(2)) , ...
                     (Estremo_2_m(1) - Estremo_1_m(1)));                    % LOCAL SDR ANGLE [rad]
        Q_m = [ cos(beta_m), sin(beta_m);
             -sin(beta_m), cos(beta_m)                                      % G2L MATRIX
            ]; 

        u_sorg = ViSorgente(Centro, Estremo_1, Estremo_2, Q', Q);           % COMPUTATION OF VELOCITY ON CENTER OF PANEL i DUE TO PANEL SOURCE j
        u_sorg_m = ViSorgente(Centro, Estremo_1_m, Estremo_2_m, Q_m', Q_m); % COMPUTATION OF VELOCITY ON CENTER OF PANEL i DUE TO MIRROR PANEL SOURCE j        
        u_vort = ViVortice(Centro, Estremo_1, Estremo_2, Q', Q);            % COMPUTATION OF VELOCITY ON CENTER OF PANEL i DUE TO PANEL VORTEX j
        u_vort_m = ViVortice(Centro, Estremo_1_m, Estremo_2_m, Q_m', Q_m);  % COMPUTATION OF VELOCITY ON CENTER OF PANEL i DUE TO MIRROR PANEL VORTEX j

        A_s_12(i,j) = (u_sorg' + u_sorg_m') * n;                            % DEFINITION OF ELEMENT i,j OF A_S SUBMATRIX
        a_v_12(i) = a_v_12(i) + (u_vort' - u_vort_m') * n;                  % DEFINITION OF ELEMENT i OF a_v 
        
         if i == 1                                                          % IF-CYCLE TO DEFINE c_s AND c_v COEFFICIENT
          c_s_1(j) = (u_sorg' + u_sorg_m') * tau1_1;
          c_v_1(j) = (u_vort' - u_vort_m') * tau1_1;
        end
        if i == N1
          c_s_N(j) = (u_sorg' + u_sorg_m') * tauN_1;
          c_v_N(j) = (u_vort' - u_vort_m') * tauN_1;
        end
    end

end
c_s_12 = c_s_1 + c_s_N;                                                     % FINAL CHARACTERIZATION OF c_s
c_v_12 = sum(c_v_1) + sum(c_v_N);                                           % FINAL CHARACTERIZATION OF c_s

% 13

c_s_1 = zeros(1,N3);
c_s_N = zeros(1,N3);
c_v_1 = zeros(1,N3);
c_v_N = zeros(1,N3);

for i = 1 : N1
    point_i = [nacapane_1.x_pane(i); nacapane_1.y_pane(i)];                 % DEFINITION OF COORDINATES OF FIRST POINT OF PANEL i
    point_ii = [nacapane_1.x_pane(i+1); nacapane_1.y_pane(i+1)];            % DEFINITION OF COORDINATES OF LAST POINT OF PANEL i
    [n, ~] = versors (point_i, point_ii);                                   % DEFINITION OF NORMAL VERSOR OF PANEL i
    center_x = (point_i(1) + point_ii(1))*0.5;                              % DEFINITION OF X-COORDINATE OF CENTER POINT OF PANEL i
    center_y = (point_i(2) + point_ii(2))*0.5;                              % DEFINITION OF X-COORDINATE OF CENTER POINT OF PANEL i
    Centro = [center_x; center_y];                                          % DEFINITION OF COORDINATES OF CENTER POINT OF PANEL i
    
    for j = 1 : N3
        Estremo_1 = [nacapane_3.x_pane(j); nacapane_3.y_pane(j)];           % DEFINITION OF COORDINATES OF FIRST POINT OF PANEL j               
        Estremo_2 = [nacapane_3.x_pane(j+1); nacapane_3.y_pane(j+1)];       % DEFINITION OF COORDINATES OF LAST POINT OF PANEL j
        beta = atan2((Estremo_2(2) - Estremo_1(2)) , ...
                     (Estremo_2(1) - Estremo_1(1)));                        % LOCAL SDR ANGLE [rad]
        Q = [ cos(beta), sin(beta);
             -sin(beta), cos(beta)                                          % G2L MATRIX
            ];                

        Estremo_1_m = [nacapane_3.x_pane_m(j); nacapane_3.y_pane_m(j)];         % DEFINITION OF COORDINATES OF FIRST POINT OF PANEL j               
        Estremo_2_m = [nacapane_3.x_pane_m(j+1); nacapane_3.y_pane_m(j+1)];     % DEFINITION OF COORDINATES OF LAST POINT OF PANEL j
        beta_m = atan2((Estremo_2_m(2) - Estremo_1_m(2)) , ...
                     (Estremo_2_m(1) - Estremo_1_m(1)));                    % LOCAL SDR ANGLE [rad]
        Q_m = [ cos(beta_m), sin(beta_m);
             -sin(beta_m), cos(beta_m)                                      % G2L MATRIX
            ]; 

        u_sorg = ViSorgente(Centro, Estremo_1, Estremo_2, Q', Q);           % COMPUTATION OF VELOCITY ON CENTER OF PANEL i DUE TO PANEL SOURCE j
        u_sorg_m = ViSorgente(Centro, Estremo_1_m, Estremo_2_m, Q_m', Q_m); % COMPUTATION OF VELOCITY ON CENTER OF PANEL i DUE TO MIRROR PANEL SOURCE j        
        u_vort = ViVortice(Centro, Estremo_1, Estremo_2, Q', Q);            % COMPUTATION OF VELOCITY ON CENTER OF PANEL i DUE TO PANEL VORTEX j
        u_vort_m = ViVortice(Centro, Estremo_1_m, Estremo_2_m, Q_m', Q_m);  % COMPUTATION OF VELOCITY ON CENTER OF PANEL i DUE TO MIRROR PANEL VORTEX j

        A_s_13(i,j) = (u_sorg' + u_sorg_m') * n;                            % DEFINITION OF ELEMENT i,j OF A_S SUBMATRIX
        a_v_13(i) = a_v_13(i) + (u_vort' - u_vort_m') * n;                  % DEFINITION OF ELEMENT i OF a_v 
        
         if i == 1                                                          % IF-CYCLE TO DEFINE c_s AND c_v COEFFICIENT
          c_s_1(j) = (u_sorg' + u_sorg_m') * tau1_1;
          c_v_1(j) = (u_vort' - u_vort_m') * tau1_1;
        end
        if i == N1
          c_s_N(j) = (u_sorg' + u_sorg_m') * tauN_1;
          c_v_N(j) = (u_vort' - u_vort_m') * tauN_1;
        end
    end

end
c_s_13 = c_s_1 + c_s_N;                                                     % FINAL CHARACTERIZATION OF c_s
c_v_13 = sum(c_v_1) + sum(c_v_N);                                           % FINAL CHARACTERIZATION OF c_s

% 21

c_s_1 = zeros(1,N1);
c_s_N = zeros(1,N1);
c_v_1 = zeros(1,N1);
c_v_N = zeros(1,N1);

for i = 1 : N2
    point_i = [nacapane_2.x_pane(i); nacapane_2.y_pane(i)];                 % DEFINITION OF COORDINATES OF FIRST POINT OF PANEL i
    point_ii = [nacapane_2.x_pane(i+1); nacapane_2.y_pane(i+1)];            % DEFINITION OF COORDINATES OF LAST POINT OF PANEL i
    [n, ~] = versors (point_i, point_ii);                                   % DEFINITION OF NORMAL VERSOR OF PANEL i
    center_x = (point_i(1) + point_ii(1))*0.5;                              % DEFINITION OF X-COORDINATE OF CENTER POINT OF PANEL i
    center_y = (point_i(2) + point_ii(2))*0.5;                              % DEFINITION OF X-COORDINATE OF CENTER POINT OF PANEL i
    Centro = [center_x; center_y];                                          % DEFINITION OF COORDINATES OF CENTER POINT OF PANEL i
    
    for j = 1 : N1
        Estremo_1 = [nacapane_1.x_pane(j); nacapane_1.y_pane(j)];           % DEFINITION OF COORDINATES OF FIRST POINT OF PANEL j               
        Estremo_2 = [nacapane_1.x_pane(j+1); nacapane_1.y_pane(j+1)];       % DEFINITION OF COORDINATES OF LAST POINT OF PANEL j
        beta = atan2((Estremo_2(2) - Estremo_1(2)) , ...
                     (Estremo_2(1) - Estremo_1(1)));                        % LOCAL SDR ANGLE [rad]
        Q = [ cos(beta), sin(beta);
             -sin(beta), cos(beta)                                          % G2L MATRIX
            ];                                                              

        Estremo_1_m = [nacapane_1.x_pane_m(j); nacapane_1.y_pane_m(j)];         % DEFINITION OF COORDINATES OF FIRST POINT OF PANEL j               
        Estremo_2_m = [nacapane_1.x_pane_m(j+1); nacapane_1.y_pane_m(j+1)];     % DEFINITION OF COORDINATES OF LAST POINT OF PANEL j
        beta_m = atan2((Estremo_2_m(2) - Estremo_1_m(2)) , ...
                     (Estremo_2_m(1) - Estremo_1_m(1)));                    % LOCAL SDR ANGLE [rad]
        Q_m = [ cos(beta_m), sin(beta_m);
             -sin(beta_m), cos(beta_m)                                      % G2L MATRIX
            ]; 

        u_sorg = ViSorgente(Centro, Estremo_1, Estremo_2, Q', Q);           % COMPUTATION OF VELOCITY ON CENTER OF PANEL i DUE TO PANEL SOURCE j
        u_sorg_m = ViSorgente(Centro, Estremo_1_m, Estremo_2_m, Q_m', Q_m); % COMPUTATION OF VELOCITY ON CENTER OF PANEL i DUE TO MIRROR PANEL SOURCE j        
        u_vort = ViVortice(Centro, Estremo_1, Estremo_2, Q', Q);            % COMPUTATION OF VELOCITY ON CENTER OF PANEL i DUE TO PANEL VORTEX j
        u_vort_m = ViVortice(Centro, Estremo_1_m, Estremo_2_m, Q_m', Q_m);  % COMPUTATION OF VELOCITY ON CENTER OF PANEL i DUE TO MIRROR PANEL VORTEX j

        A_s_21(i,j) = (u_sorg' + u_sorg_m') * n;                            % DEFINITION OF ELEMENT i,j OF A_S SUBMATRIX
        a_v_21(i) = a_v_21(i) + (u_vort' - u_vort_m') * n;                  % DEFINITION OF ELEMENT i OF a_v 
        
        if i == 1                                                           % IF-CYCLE TO DEFINE c_s AND c_v COEFFICIENT
          c_s_1(j) = (u_sorg' + u_sorg_m') * tau1_2;
          c_v_1(j) = (u_vort' - u_vort_m') * tau1_2;
        end
        if i == N2
          c_s_N(j) = (u_sorg' + u_sorg_m') * tauN_2;
          c_v_N(j) = (u_vort' - u_vort_m') * tauN_2;
        end
    end

end
c_s_21 = c_s_1 + c_s_N;                                                     % FINAL CHARACTERIZATION OF c_s
c_v_21 = sum(c_v_1) + sum(c_v_N);                                           % FINAL CHARACTERIZATION OF c_s

% 22

c_s_1 = zeros(1,N2);
c_s_N = zeros(1,N2);
c_v_1 = zeros(1,N2);
c_v_N = zeros(1,N2);

for i = 1 : N2
    point_i = [nacapane_2.x_pane(i); nacapane_2.y_pane(i)];                 % DEFINITION OF COORDINATES OF FIRST POINT OF PANEL i
    point_ii = [nacapane_2.x_pane(i+1); nacapane_2.y_pane(i+1)];            % DEFINITION OF COORDINATES OF LAST POINT OF PANEL i
    [n, ~] = versors (point_i, point_ii);                                   % DEFINITION OF NORMAL VERSOR OF PANEL i
    center_x = (point_i(1) + point_ii(1))*0.5;                              % DEFINITION OF X-COORDINATE OF CENTER POINT OF PANEL i
    center_y = (point_i(2) + point_ii(2))*0.5;                              % DEFINITION OF X-COORDINATE OF CENTER POINT OF PANEL i
    Centro = [center_x; center_y];                                          % DEFINITION OF COORDINATES OF CENTER POINT OF PANEL i
    center_2(i, :) = Centro;                                                % STORING THE COORDINATE OF CENTER POINTS OF ALL PANELS

    for j = 1 : N2
        Estremo_1 = [nacapane_2.x_pane(j); nacapane_2.y_pane(j)];           % DEFINITION OF COORDINATES OF FIRST POINT OF PANEL j               
        Estremo_2 = [nacapane_2.x_pane(j+1); nacapane_2.y_pane(j+1)];       % DEFINITION OF COORDINATES OF LAST POINT OF PANEL j
        beta = atan2((Estremo_2(2) - Estremo_1(2)) , ...
                     (Estremo_2(1) - Estremo_1(1)));                        % LOCAL SDR ANGLE [rad]
        Q = [ cos(beta), sin(beta);
             -sin(beta), cos(beta)                                          % G2L MATRIX
            ];                                                              

        Estremo_1_m = [nacapane_2.x_pane_m(j); nacapane_2.y_pane_m(j)];         % DEFINITION OF COORDINATES OF FIRST POINT OF PANEL j               
        Estremo_2_m = [nacapane_2.x_pane_m(j+1); nacapane_2.y_pane_m(j+1)];     % DEFINITION OF COORDINATES OF LAST POINT OF PANEL j
        beta_m = atan2((Estremo_2_m(2) - Estremo_1_m(2)) , ...
                     (Estremo_2_m(1) - Estremo_1_m(1)));                    % LOCAL SDR ANGLE [rad]
        Q_m = [ cos(beta_m), sin(beta_m);
             -sin(beta_m), cos(beta_m)                                      % G2L MATRIX
            ]; 

        u_sorg = ViSorgente(Centro, Estremo_1, Estremo_2, Q', Q);           % COMPUTATION OF VELOCITY ON CENTER OF PANEL i DUE TO PANEL SOURCE j
        u_sorg_m = ViSorgente(Centro, Estremo_1_m, Estremo_2_m, Q_m', Q_m); % COMPUTATION OF VELOCITY ON CENTER OF PANEL i DUE TO MIRROR PANEL SOURCE j        
        u_vort = ViVortice(Centro, Estremo_1, Estremo_2, Q', Q);            % COMPUTATION OF VELOCITY ON CENTER OF PANEL i DUE TO PANEL VORTEX j
        u_vort_m = ViVortice(Centro, Estremo_1_m, Estremo_2_m, Q_m', Q_m);  % COMPUTATION OF VELOCITY ON CENTER OF PANEL i DUE TO MIRROR PANEL VORTEX j

        A_s_22(i,j) = (u_sorg' + u_sorg_m') * n;                            % DEFINITION OF ELEMENT i,j OF A_S SUBMATRIX
        a_v_22(i) = a_v_22(i) + (u_vort' - u_vort_m') * n;                  % DEFINITION OF ELEMENT i OF a_v 
        
         if i == 1                                                          % IF-CYCLE TO DEFINE c_s AND c_v COEFFICIENT
          c_s_1(j) = (u_sorg' + u_sorg_m') * tau1_2;
          c_v_1(j) = (u_vort' - u_vort_m') * tau1_2;
        end
        if i == N2
          c_s_N(j) = (u_sorg' + u_sorg_m') * tauN_2;
          c_v_N(j) = (u_vort' - u_vort_m') * tauN_2;
        end

    end

    b_s_2(i) = -u_inf' * n;                                                 % DEFINITION OF b_s VECTOR
    b_v_2 = - u_inf' * (tau1_2 + tauN_2);                                   % DEFINITION OF b_v VARIABLE

end
c_s_22 = c_s_1 + c_s_N;                                                     % FINAL CHARACTERIZATION OF c_s
c_v_22 = sum(c_v_1) + sum(c_v_N);                                           % FINAL CHARACTERIZATION OF c_s

% 23

c_s_1 = zeros(1,N3);
c_s_N = zeros(1,N3);
c_v_1 = zeros(1,N3);
c_v_N = zeros(1,N3);

for i = 1 : N2
    point_i = [nacapane_2.x_pane(i); nacapane_2.y_pane(i)];                 % DEFINITION OF COORDINATES OF FIRST POINT OF PANEL i
    point_ii = [nacapane_2.x_pane(i+1); nacapane_2.y_pane(i+1)];            % DEFINITION OF COORDINATES OF LAST POINT OF PANEL i
    [n, ~] = versors (point_i, point_ii);                                   % DEFINITION OF NORMAL VERSOR OF PANEL i
    center_x = (point_i(1) + point_ii(1))*0.5;                              % DEFINITION OF X-COORDINATE OF CENTER POINT OF PANEL i
    center_y = (point_i(2) + point_ii(2))*0.5;                              % DEFINITION OF X-COORDINATE OF CENTER POINT OF PANEL i
    Centro = [center_x; center_y];                                          % DEFINITION OF COORDINATES OF CENTER POINT OF PANEL i
    
    for j = 1 : N3
        Estremo_1 = [nacapane_3.x_pane(j); nacapane_3.y_pane(j)];           % DEFINITION OF COORDINATES OF FIRST POINT OF PANEL j               
        Estremo_2 = [nacapane_3.x_pane(j+1); nacapane_3.y_pane(j+1)];       % DEFINITION OF COORDINATES OF LAST POINT OF PANEL j
        beta = atan2((Estremo_2(2) - Estremo_1(2)) , ...
                     (Estremo_2(1) - Estremo_1(1)));                        % LOCAL SDR ANGLE [rad]
        Q = [ cos(beta), sin(beta);
             -sin(beta), cos(beta)                                          % G2L MATRIX
            ];                                                              

        Estremo_1_m = [nacapane_3.x_pane_m(j); nacapane_3.y_pane_m(j)];         % DEFINITION OF COORDINATES OF FIRST POINT OF PANEL j               
        Estremo_2_m = [nacapane_3.x_pane_m(j+1); nacapane_3.y_pane_m(j+1)];     % DEFINITION OF COORDINATES OF LAST POINT OF PANEL j
        beta_m = atan2((Estremo_2_m(2) - Estremo_1_m(2)) , ...
                     (Estremo_2_m(1) - Estremo_1_m(1)));                    % LOCAL SDR ANGLE [rad]
        Q_m = [ cos(beta_m), sin(beta_m);
             -sin(beta_m), cos(beta_m)                                      % G2L MATRIX
            ]; 

        u_sorg = ViSorgente(Centro, Estremo_1, Estremo_2, Q', Q);           % COMPUTATION OF VELOCITY ON CENTER OF PANEL i DUE TO PANEL SOURCE j
        u_sorg_m = ViSorgente(Centro, Estremo_1_m, Estremo_2_m, Q_m', Q_m); % COMPUTATION OF VELOCITY ON CENTER OF PANEL i DUE TO MIRROR PANEL SOURCE j        
        u_vort = ViVortice(Centro, Estremo_1, Estremo_2, Q', Q);            % COMPUTATION OF VELOCITY ON CENTER OF PANEL i DUE TO PANEL VORTEX j
        u_vort_m = ViVortice(Centro, Estremo_1_m, Estremo_2_m, Q_m', Q_m);  % COMPUTATION OF VELOCITY ON CENTER OF PANEL i DUE TO MIRROR PANEL VORTEX j

        A_s_23(i,j) = (u_sorg' + u_sorg_m') * n;                            % DEFINITION OF ELEMENT i,j OF A_S SUBMATRIX
        a_v_23(i) = a_v_23(i) + (u_vort' - u_vort_m') * n;                  % DEFINITION OF ELEMENT i OF a_v 
        
        if i == 1                                                           % IF-CYCLE TO DEFINE c_s AND c_v COEFFICIENT
          c_s_1(j) = (u_sorg' + u_sorg_m') * tau1_2;
          c_v_1(j) = (u_vort' - u_vort_m') * tau1_2;
        end
        if i == N2
          c_s_N(j) = (u_sorg' + u_sorg_m') * tauN_2;
          c_v_N(j) = (u_vort' - u_vort_m') * tauN_2;
        end
    end

end
c_s_23 = c_s_1 + c_s_N;                                                     % FINAL CHARACTERIZATION OF c_s
c_v_23 = sum(c_v_1) + sum(c_v_N);                                           % FINAL CHARACTERIZATION OF c_s

% 31

c_s_1 = zeros(1,N1);
c_s_N = zeros(1,N1);
c_v_1 = zeros(1,N1);
c_v_N = zeros(1,N1);

for i = 1 : N3
    point_i = [nacapane_3.x_pane(i); nacapane_3.y_pane(i)];                 % DEFINITION OF COORDINATES OF FIRST POINT OF PANEL i
    point_ii = [nacapane_3.x_pane(i+1); nacapane_3.y_pane(i+1)];            % DEFINITION OF COORDINATES OF LAST POINT OF PANEL i
    [n, ~] = versors (point_i, point_ii);                                   % DEFINITION OF NORMAL VERSOR OF PANEL i
    center_x = (point_i(1) + point_ii(1))*0.5;                              % DEFINITION OF X-COORDINATE OF CENTER POINT OF PANEL i
    center_y = (point_i(2) + point_ii(2))*0.5;                              % DEFINITION OF X-COORDINATE OF CENTER POINT OF PANEL i
    Centro = [center_x; center_y];                                          % DEFINITION OF COORDINATES OF CENTER POINT OF PANEL i
    
    for j = 1 : N1
        Estremo_1 = [nacapane_1.x_pane(j); nacapane_1.y_pane(j)];           % DEFINITION OF COORDINATES OF FIRST POINT OF PANEL j               
        Estremo_2 = [nacapane_1.x_pane(j+1); nacapane_1.y_pane(j+1)];       % DEFINITION OF COORDINATES OF LAST POINT OF PANEL j
        beta = atan2((Estremo_2(2) - Estremo_1(2)) , ...
                     (Estremo_2(1) - Estremo_1(1)));                        % LOCAL SDR ANGLE [rad]
        Q = [ cos(beta), sin(beta);
             -sin(beta), cos(beta)                                          % G2L MATRIX
            ];                                                              

        Estremo_1_m = [nacapane_1.x_pane_m(j); nacapane_1.y_pane_m(j)];         % DEFINITION OF COORDINATES OF FIRST POINT OF PANEL j               
        Estremo_2_m = [nacapane_1.x_pane_m(j+1); nacapane_1.y_pane_m(j+1)];     % DEFINITION OF COORDINATES OF LAST POINT OF PANEL j
        beta_m = atan2((Estremo_2_m(2) - Estremo_1_m(2)) , ...
                     (Estremo_2_m(1) - Estremo_1_m(1)));                    % LOCAL SDR ANGLE [rad]
        Q_m = [ cos(beta_m), sin(beta_m);
             -sin(beta_m), cos(beta_m)                                      % G2L MATRIX
            ]; 

        u_sorg = ViSorgente(Centro, Estremo_1, Estremo_2, Q', Q);           % COMPUTATION OF VELOCITY ON CENTER OF PANEL i DUE TO PANEL SOURCE j
        u_sorg_m = ViSorgente(Centro, Estremo_1_m, Estremo_2_m, Q_m', Q_m); % COMPUTATION OF VELOCITY ON CENTER OF PANEL i DUE TO MIRROR PANEL SOURCE j        
        u_vort = ViVortice(Centro, Estremo_1, Estremo_2, Q', Q);            % COMPUTATION OF VELOCITY ON CENTER OF PANEL i DUE TO PANEL VORTEX j
        u_vort_m = ViVortice(Centro, Estremo_1_m, Estremo_2_m, Q_m', Q_m);  % COMPUTATION OF VELOCITY ON CENTER OF PANEL i DUE TO MIRROR PANEL VORTEX j

        A_s_31(i,j) = (u_sorg' + u_sorg_m') * n;                            % DEFINITION OF ELEMENT i,j OF A_S SUBMATRIX
        a_v_31(i) = a_v_31(i) + (u_vort' - u_vort_m') * n;                  % DEFINITION OF ELEMENT i OF a_v 
        
        if i == 1                                                           % IF-CYCLE TO DEFINE c_s AND c_v COEFFICIENT
          c_s_1(j) = (u_sorg' + u_sorg_m') * tau1_3;
          c_v_1(j) = (u_vort' - u_vort_m') * tau1_3;
        end
        if i == N3
          c_s_N(j) = (u_sorg' + u_sorg_m') * tauN_3;
          c_v_N(j) = (u_vort' - u_vort_m') * tauN_3;
        end
    end

end
c_s_31 = c_s_1 + c_s_N;                                                     % FINAL CHARACTERIZATION OF c_s
c_v_31 = sum(c_v_1) + sum(c_v_N);                                           % FINAL CHARACTERIZATION OF c_s

% 32

c_s_1 = zeros(1,N2);
c_s_N = zeros(1,N2);
c_v_1 = zeros(1,N2);
c_v_N = zeros(1,N2);

for i = 1 : N3
    point_i = [nacapane_3.x_pane(i); nacapane_3.y_pane(i)];                 % DEFINITION OF COORDINATES OF FIRST POINT OF PANEL i
    point_ii = [nacapane_3.x_pane(i+1); nacapane_3.y_pane(i+1)];            % DEFINITION OF COORDINATES OF LAST POINT OF PANEL i
    [n, ~] = versors (point_i, point_ii);                                   % DEFINITION OF NORMAL VERSOR OF PANEL i
    center_x = (point_i(1) + point_ii(1))*0.5;                              % DEFINITION OF X-COORDINATE OF CENTER POINT OF PANEL i
    center_y = (point_i(2) + point_ii(2))*0.5;                              % DEFINITION OF X-COORDINATE OF CENTER POINT OF PANEL i
    Centro = [center_x; center_y];                                          % DEFINITION OF COORDINATES OF CENTER POINT OF PANEL i
    
    for j = 1 : N2
        Estremo_1 = [nacapane_2.x_pane(j); nacapane_2.y_pane(j)];           % DEFINITION OF COORDINATES OF FIRST POINT OF PANEL j               
        Estremo_2 = [nacapane_2.x_pane(j+1); nacapane_2.y_pane(j+1)];       % DEFINITION OF COORDINATES OF LAST POINT OF PANEL j
        beta = atan2((Estremo_2(2) - Estremo_1(2)) , ...
                     (Estremo_2(1) - Estremo_1(1)));                        % LOCAL SDR ANGLE [rad]
        Q = [ cos(beta), sin(beta);
             -sin(beta), cos(beta)                                          % G2L MATRIX
            ];                                                              

        Estremo_1_m = [nacapane_2.x_pane_m(j); nacapane_2.y_pane_m(j)];         % DEFINITION OF COORDINATES OF FIRST POINT OF PANEL j               
        Estremo_2_m = [nacapane_2.x_pane_m(j+1); nacapane_2.y_pane_m(j+1)];     % DEFINITION OF COORDINATES OF LAST POINT OF PANEL j
        beta_m = atan2((Estremo_2_m(2) - Estremo_1_m(2)) , ...
                     (Estremo_2_m(1) - Estremo_1_m(1)));                    % LOCAL SDR ANGLE [rad]
        Q_m = [ cos(beta_m), sin(beta_m);
             -sin(beta_m), cos(beta_m)                                      % G2L MATRIX
            ]; 

        u_sorg = ViSorgente(Centro, Estremo_1, Estremo_2, Q', Q);           % COMPUTATION OF VELOCITY ON CENTER OF PANEL i DUE TO PANEL SOURCE j
        u_sorg_m = ViSorgente(Centro, Estremo_1_m, Estremo_2_m, Q_m', Q_m); % COMPUTATION OF VELOCITY ON CENTER OF PANEL i DUE TO MIRROR PANEL SOURCE j        
        u_vort = ViVortice(Centro, Estremo_1, Estremo_2, Q', Q);            % COMPUTATION OF VELOCITY ON CENTER OF PANEL i DUE TO PANEL VORTEX j
        u_vort_m = ViVortice(Centro, Estremo_1_m, Estremo_2_m, Q_m', Q_m);  % COMPUTATION OF VELOCITY ON CENTER OF PANEL i DUE TO MIRROR PANEL VORTEX j

        A_s_32(i,j) = (u_sorg' + u_sorg_m') * n;                            % DEFINITION OF ELEMENT i,j OF A_S SUBMATRIX
        a_v_32(i) = a_v_32(i) + (u_vort' - u_vort_m') * n;                  % DEFINITION OF ELEMENT i OF a_v 
        
        if i == 1                                                           % IF-CYCLE TO DEFINE c_s AND c_v COEFFICIENT
          c_s_1(j) = (u_sorg' + u_sorg_m') * tau1_3;
          c_v_1(j) = (u_vort' - u_vort_m') * tau1_3;
        end
        if i == N3
          c_s_N(j) = (u_sorg' + u_sorg_m') * tauN_3;
          c_v_N(j) = (u_vort' - u_vort_m') * tauN_3;
        end
    end

end
c_s_32 = c_s_1 + c_s_N;                                                     % FINAL CHARACTERIZATION OF c_s
c_v_32 = sum(c_v_1) + sum(c_v_N);                                           % FINAL CHARACTERIZATION OF c_s

% 33

c_s_1 = zeros(1,N3);
c_s_N = zeros(1,N3);
c_v_1 = zeros(1,N3);
c_v_N = zeros(1,N3);

for i = 1 : N3
    point_i = [nacapane_3.x_pane(i); nacapane_3.y_pane(i)];                 % DEFINITION OF COORDINATES OF FIRST POINT OF PANEL i
    point_ii = [nacapane_3.x_pane(i+1); nacapane_3.y_pane(i+1)];            % DEFINITION OF COORDINATES OF LAST POINT OF PANEL i
    [n, ~] = versors (point_i, point_ii);                                   % DEFINITION OF NORMAL VERSOR OF PANEL i
    center_x = (point_i(1) + point_ii(1))*0.5;                              % DEFINITION OF X-COORDINATE OF CENTER POINT OF PANEL i
    center_y = (point_i(2) + point_ii(2))*0.5;                              % DEFINITION OF X-COORDINATE OF CENTER POINT OF PANEL i
    Centro = [center_x; center_y];                                          % DEFINITION OF COORDINATES OF CENTER POINT OF PANEL i
    center_3(i, :) = Centro;                                                % STORING THE COORDINATE OF CENTER POINTS OF ALL PANELS

    for j = 1 : N3
        Estremo_1 = [nacapane_3.x_pane(j); nacapane_3.y_pane(j)];           % DEFINITION OF COORDINATES OF FIRST POINT OF PANEL j               
        Estremo_2 = [nacapane_3.x_pane(j+1); nacapane_3.y_pane(j+1)];       % DEFINITION OF COORDINATES OF LAST POINT OF PANEL j
        beta = atan2((Estremo_2(2) - Estremo_1(2)) , ...
                     (Estremo_2(1) - Estremo_1(1)));                        % LOCAL SDR ANGLE [rad]
        Q = [ cos(beta), sin(beta);
             -sin(beta), cos(beta)                                          % G2L MATRIX
            ];                                                              

        Estremo_1_m = [nacapane_3.x_pane_m(j); nacapane_3.y_pane_m(j)];         % DEFINITION OF COORDINATES OF FIRST POINT OF PANEL j               
        Estremo_2_m = [nacapane_3.x_pane_m(j+1); nacapane_3.y_pane_m(j+1)];     % DEFINITION OF COORDINATES OF LAST POINT OF PANEL j
        beta_m = atan2((Estremo_2_m(2) - Estremo_1_m(2)) , ...
                     (Estremo_2_m(1) - Estremo_1_m(1)));                    % LOCAL SDR ANGLE [rad]
        Q_m = [ cos(beta_m), sin(beta_m);
             -sin(beta_m), cos(beta_m)                                      % G2L MATRIX
            ]; 

        u_sorg = ViSorgente(Centro, Estremo_1, Estremo_2, Q', Q);           % COMPUTATION OF VELOCITY ON CENTER OF PANEL i DUE TO PANEL SOURCE j
        u_sorg_m = ViSorgente(Centro, Estremo_1_m, Estremo_2_m, Q_m', Q_m); % COMPUTATION OF VELOCITY ON CENTER OF PANEL i DUE TO MIRROR PANEL SOURCE j        
        u_vort = ViVortice(Centro, Estremo_1, Estremo_2, Q', Q);            % COMPUTATION OF VELOCITY ON CENTER OF PANEL i DUE TO PANEL VORTEX j
        u_vort_m = ViVortice(Centro, Estremo_1_m, Estremo_2_m, Q_m', Q_m);  % COMPUTATION OF VELOCITY ON CENTER OF PANEL i DUE TO MIRROR PANEL VORTEX j

        A_s_33(i,j) = (u_sorg' + u_sorg_m') * n;                            % DEFINITION OF ELEMENT i,j OF A_S SUBMATRIX
        a_v_33(i) = a_v_33(i) + (u_vort' - u_vort_m') * n;                  % DEFINITION OF ELEMENT i OF a_v 
        
         if i == 1                                                          % IF-CYCLE TO DEFINE c_s AND c_v COEFFICIENT
          c_s_1(j) = (u_sorg' + u_sorg_m') * tau1_3;
          c_v_1(j) = (u_vort' - u_vort_m') * tau1_3;
        end
        if i == N3
          c_s_N(j) = (u_sorg' + u_sorg_m') * tauN_3;
          c_v_N(j) = (u_vort' - u_vort_m') * tauN_3;
        end

    end

    b_s_3(i) = -u_inf' * n;                                                 % DEFINITION OF b_s VECTOR
    b_v_3 = - u_inf' * (tau1_3 + tauN_3);                                   % DEFINITION OF b_v VARIABLE

end
c_s_33 = c_s_1 + c_s_N;                                                     % FINAL CHARACTERIZATION OF c_s
c_v_33 = sum(c_v_1) + sum(c_v_N);                                           % FINAL CHARACTERIZATION OF c_s

%% MATRIX DEFINITION
A = [A_s_11, A_s_12, A_s_13, a_v_11, a_v_12, a_v_13;
     A_s_21, A_s_22, A_s_23, a_v_21, a_v_22, a_v_23;
     A_s_31, A_s_32, A_s_33, a_v_31, a_v_32, a_v_33;
     c_s_11, c_s_12, c_s_13, c_v_11, c_v_12, c_v_13;
     c_s_21, c_s_22, c_s_23, c_v_21, c_v_22, c_v_23;
     c_s_31, c_s_32, c_s_33, c_v_31, c_v_32, c_v_33];                       % ASSEMBLY OF LINEAR SYSTEM MATRIX

b = [b_s_1;
     b_s_2;
     b_s_3;
     b_v_1;
     b_v_2;
     b_v_3];                                                                % ASSEMBLY OF CONSTANT TERM