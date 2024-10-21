%% INTRO
clc;
clear;
close all;

%% PROBLEM DEFINITION

run("profili.m")

nacapane_1 = nacapane_slat;
N1 = N_slat;

nacapane_2 = nacapane_main;
N2 = N_main;

nacapane_3 = nacapane_flap;
N3 = N_flap;

u_inf = [50;
          0];                                                               % ASINTOTIC VELOCITY [m/s]


%% FOIL PANELISATION PLOT 
% figure;
axis equal;
hold on;
grid on;
plot(nacapane_1.x_pane, nacapane_1.y_pane, '-');                            % FOIL 1 RAPRESENTATION
plot(nacapane_2.x_pane, nacapane_2.y_pane, '-');                            % FOIL 2 RAPRESENTATION
plot(nacapane_3.x_pane, nacapane_3.y_pane, '-');                            % FOIL 2 RAPRESENTATION

plot(nacapane_1.x_pane_m, nacapane_1.y_pane_m, '-');                        % FOIL 1 RAPRESENTATION
plot(nacapane_2.x_pane_m, nacapane_2.y_pane_m, '-');                        % FOIL 2 RAPRESENTATION
plot(nacapane_3.x_pane_m, nacapane_3.y_pane_m, '-');                        % FOIL 2 RAPRESENTATION

% plot(center_1(:, 1), center_1(:, 2), 'o');                                % CENTERS 1
% plot(center_2(:, 1), center_2(:, 2), 'o');                                % CENTERS 2
% plot(center_3(:, 1), center_3(:, 2), 'o');                                % CENTERS 2

title("Foil panelisation", 'Interpreter', 'latex')
legend('slat', 'main', 'flap')

%% LINEAR SISTEM
[A, b, center_1, center_2, center_3] = Matrix_definition(nacapane_1, nacapane_2, nacapane_3, N1, N2, N3, u_inf);

x = A\b;                                                                    % COMPUTATION OF LINEAR SISTEM SOLUTION
q_1 = x(1:N1);                                                              % INTENSITIES OF SURCES FOIL 1
q_2 = x(N1+1:N1+N2);                                                        % INTENSITIES OF SURCES FOIL 2
q_3 = x(N1+N2+1:end-3);                                                        % INTENSITIES OF SURCES FOIL 3

gamma_1 = x(end-2);                                                         % INTENSITY OF VORTEX FOIL 1
gamma_2 = x(end-1);                                                         % INTENSITY OF VORTEX FOIL 2
gamma_3 = x(end);                                                           % INTENSITY OF VORTEX FOIL 3

u_tot = @(x, y) TotalVelocity_groundX3 (q_1, q_2, q_3, gamma_1, gamma_2, ...
                  gamma_3, u_inf, x, y, nacapane_1, nacapane_2, ...
                  nacapane_3, N1, N2, N3);


%% plot y-component of velocity
Nfx = 100; 
Nfy = 100;

xx = linspace(-1,4,Nfx);
V_y_field = zeros(1, Nfy);

for i = 1:Nfx
        u_aux = TotalVelocity_groundX3 (q_1, q_2, q_3,gamma_1, gamma_2, ...
                                        gamma_3, u_inf, xx(i), 0, ...
                                        nacapane_1, nacapane_2, ...
                                        nacapane_3, N1, N2, N3);
        V_y_field(i) = u_aux'*[0;1];
end

figure;
plot(xx, V_y_field);


%% COEFFICIENTS COMPUTATION

Cp_1 = zeros(N1, 1);                                                        % MEMORY ALLOCATION
Cp_2 = zeros(N2, 1);                                                        % MEMORY ALLOCATION
Cp_3 = zeros(N3, 1);                                                        % MEMORY ALLOCATION

Circ_1 = 0;                                                                 % VARIABLE INITIALIZATION
Circ_2 = 0;                                                                 % VARIABLE INITIALIZATION
Circ_3 = 0;                                                                 % VARIABLE INITIALIZATION

% 1
for i = 1 : N1

    Estremo_1 = [nacapane_1.x_pane(i); nacapane_1.y_pane(i)];               % DEFINITION OF COORDINATES OF FIRST POINT OF PANEL i           
    Estremo_2 = [nacapane_1.x_pane(i+1); nacapane_1.y_pane(i+1)];           % DEFINITION OF COORDINATES OF LAST POINT OF PANEL i           

    [~, tau] = versors (Estremo_1, Estremo_2);                              % DEFINITION OF TANGENT VERSOR OF PANEL i
    
    aux = u_tot(center_1(i,1),center_1(i,2))'*tau;
    Cp_1(i) = 1 - aux^2/(norm(u_inf))^2;                                    % PRESSURE COEFFICIENT COMPUTATION

    Circ_1 = Circ_1 + norm(Estremo_2 - Estremo_1)*gamma_1;                  % CIRCOLATION COMPUTATION
end

Cl_1 = 2*Circ_1/(norm(u_inf)*chord_slat);                                                % LIFT COEFFICIENT COMPUTATION

% 2

for i = 1 : N2

    Estremo_1 = [nacapane_2.x_pane(i); nacapane_2.y_pane(i)];               % DEFINITION OF COORDINATES OF FIRST POINT OF PANEL i           
    Estremo_2 = [nacapane_2.x_pane(i+1); nacapane_2.y_pane(i+1)];           % DEFINITION OF COORDINATES OF LAST POINT OF PANEL i           

    [~, tau] = versors (Estremo_1, Estremo_2);                              % DEFINITION OF TANGENT VERSOR OF PANEL i
    
    aux = u_tot(center_2(i,1),center_2(i,2))'*tau;
    Cp_2(i) = 1 - aux^2/(norm(u_inf))^2;                                    % PRESSURE COEFFICIENT COMPUTATION

    Circ_2 = Circ_2 + norm(Estremo_2 - Estremo_1)*gamma_2;                  % CIRCOLATION COMPUTATION
end

Cl_2 = 2*Circ_2/(norm(u_inf)*chord_main);                                                % LIFT COEFFICIENT COMPUTATION

% 3
for i = 1 : N3

    Estremo_1 = [nacapane_3.x_pane(i); nacapane_3.y_pane(i)];               % DEFINITION OF COORDINATES OF FIRST POINT OF PANEL i           
    Estremo_2 = [nacapane_3.x_pane(i+1); nacapane_3.y_pane(i+1)];           % DEFINITION OF COORDINATES OF LAST POINT OF PANEL i           

    [~, tau] = versors (Estremo_1, Estremo_2);                              % DEFINITION OF TANGENT VERSOR OF PANEL i
    
    aux = u_tot(center_3(i,1),center_3(i,2))'*tau;
    Cp_3(i) = 1 - aux^2/(norm(u_inf))^2;                                    % PRESSURE COEFFICIENT COMPUTATION

    Circ_3 = Circ_3 + norm(Estremo_2 - Estremo_1)*gamma_3;                  % CIRCOLATION COMPUTATION
end

Cl_3 = 2*Circ_3/(norm(u_inf)*chord_flap);                                                % LIFT COEFFICIENT COMPUTATION

rho = 1.225;
L = rho*norm(u_inf)*(Circ_1+Circ_2+Circ_3);
Cl = L/(0.5*rho*norm(u_inf)^2*(chord_main+chord_slat+chord_flap));
%% PRESSURE COEFFICIENT PLOT
figure;
hold on;
grid on;
xplot = linspace(0,1,N1/2+1);
xplot_meno = linspace(1,0,N1/2);
plot(xplot_meno,-Cp_1(1:N1/2), 'b-', 'LineWidth', 1.5)
plot(xplot,-Cp_1(N1/2:N1), 'b--', 'LineWidth', 1.5)
plot(xplot_meno,-Cp_2(1:N1/2), 'r-', 'LineWidth', 1.5)     
plot(xplot,-Cp_2(N2/2:N2), 'r--', 'LineWidth', 1.5)

plot(xplot_meno,-Cp_3(1:N3/2), 'g-', 'LineWidth', 1.5)     
plot(xplot,-Cp_3(N3/2:N3), 'g--', 'LineWidth', 1.5)

legend('bottom airfoil 1', 'top airfoil 1', ...
       'bottom airfoil 2', 'top airfoil 2', ...
       'bottom airfoil 3', 'top airfoil 3')
xlabel('$x$', 'Interpreter', 'latex')
ylabel('$C_p$','Interpreter', 'latex')

figure
hold on
grid on
plot(nacapane_1.x_pane(1:N1/2),-Cp_1(1:N1/2), 'b-', 'LineWidth', 1.5)
plot(nacapane_1.x_pane(N1/2:N1),-Cp_1(N1/2:N1), 'b--', 'LineWidth', 1.5)
plot(nacapane_2.x_pane(1:N2/2) ,-Cp_2(1:N2/2), 'r-', 'LineWidth', 1.5)     
plot(nacapane_2.x_pane(N2/2:N2) ,-Cp_2(N2/2:N2), 'r--', 'LineWidth', 1.5)
plot(nacapane_3.x_pane(1:N3/2) ,-Cp_3(1:N3/2), 'g-', 'LineWidth', 1.5)     
plot(nacapane_3.x_pane(N3/2:N3) ,-Cp_3(N3/2:N3), 'g--', 'LineWidth', 1.5)



%% FIELD PLOT
x_lim = [-1, 2];
y_lim = [0, h_main+1];
Nfx = 50;
Nfy = 50;

Velocity_Plot(Nfx, Nfy, x_lim, y_lim, q_1, q_2, q_3, gamma_1, gamma_2, ...
              gamma_3, nacapane_1, nacapane_2, nacapane_3, N1, N2, N3, ...
              u_inf, 'A_Field');
 