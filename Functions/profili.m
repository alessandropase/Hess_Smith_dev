%% INTRO

close all;
clc;
clear;

%% PROBLEM DEFINITION

FileName_1 = 'naca_0012_closed';                                            % FILE NAME MAIN FOIL
FileName_2 = 'naca_23012_closed';                                           % FILE NAME FLAP
FileName_3 = 'naca_9315_closed';                                            % FILE NAME SLAT

chord_main = 1;                                                             % Chord length of the main airfoil
chord_flap = 0.20 * chord_main;                                             % Chord length of the flap, expressed in percentage of main chord
chord_slat = 0.15* chord_main;                                              % Chord length of the slat, expressed in percentage of main chord

h_main = 0.7;                                                                 % Distance from the ground [m]

%% TAKE-OFF CONFIG

% alpha_main_deg =  8;                                                        % AoA [deg]
% alpha_main = alpha_main_deg * pi / 180;                                     % AoA [rad]
% alpha_deg_flap = 20 + alpha_main_deg;                                       % AoA [deg]
% alpha_flap = alpha_deg_flap * pi / 180;                                     % AoA [rad]
% alpha_deg_slat = -35;                                                       % AoA [deg]
% alpha_slat = alpha_deg_slat * pi / 180;                                     % AoA [rad]
% 
% gap_x_flap2main = 0.005;
% gap_y_flap2main = -0.05;
% gap_x_slat2main = -0.35;
% gap_y_slat2main = 0.35;


%% LANDING CONFIG

alpha_main_deg =  0;                                                        % AoA [deg]
alpha_main = alpha_main_deg * pi / 180;                                     % AoA [rad]
alpha_deg_flap = 30;                                                        % AoA [deg]
alpha_flap = alpha_deg_flap * pi / 180;                                     % AoA [rad]
alpha_deg_slat = -30;                                                       % AoA [deg]
alpha_slat = alpha_deg_slat * pi / 180;                                     % AoA [rad]

gap_x_flap2main = 0.015;
gap_y_flap2main = -0.025;
gap_x_slat2main = -0.09;
gap_y_slat2main = -0.008;

% gap_x_flap2main = 5;
% gap_y_flap2main = 5;
% gap_x_slat2main = 5;
% gap_y_slat2main = 5;

%% PANELISATION OF FOILS
% MAIN
nacapane_tab_main = importPanelisation(FileName_1);                         % IMPORT OF PANELISATION FROM FILE .txt
array_main = table2array(nacapane_tab_main);                                % CONVERTION OF THE TABLE TO AN ARRAY

array_main = chord_main*ROT_alpha(array_main, alpha_main);                  % ROTATION OF THE FOIL OF AoA

h_main = h_main + abs(min(array_main(:, 2)));
array_main = flip(array_main);

nacapane_main.x_pane = array_main(:, 1);                                    % DEFINING THE STRUCT OF POINTS - X
nacapane_main.y_pane = array_main(:, 2) + h_main;                           % DEFINING THE STRUCT OF POINTS - Y

N_main = length(nacapane_main.x_pane) - 1;                                  % NUMBER OF PANELS

% FLAP
nacapane_tab_flap = importPanelisation(FileName_2);                         % IMPORT OF PANELISATION FROM FILE .txt
array_flap = table2array(nacapane_tab_flap);                                % CONVERTION OF THE TABLE TO AN ARRAY
                                      
array_flap = ROT_alpha(array_flap, alpha_flap);                             % ROTATION OF THE FOIL OF AoA
% In questo caso non tolgo 0.25 perchè lo voglio ruotare attorno al suo l.e.
array_flap = chord_flap*flip(array_flap);
offset_x_flap = gap_x_flap2main + nacapane_main.x_pane(1);                  % FOIL 2 X-OFFSET
offset_y_flap = gap_y_flap2main + nacapane_main.y_pane(1);
nacapane_flap.x_pane = array_flap(:, 1) + offset_x_flap;                    % DEFINING THE STRUCT OF POINTS - X
nacapane_flap.y_pane = array_flap(:, 2) + offset_y_flap;                    % DEFINING THE STRUCT OF POINTS - Y

N_flap = length(nacapane_flap.x_pane) - 1;                                  % NUMBER OF PANELS

% SLAT 
nacapane_tab_slat = importPanelisation(FileName_3);                         % IMPORT OF PANELISATION FROM FILE .txt
array_slat = table2array(nacapane_tab_slat);                                % CONVERTION OF THE TABLE TO AN ARRAY

array_slat = ROT_alpha(array_slat, alpha_slat);                             % ROTATION OF THE FOIL OF AoA
array_slat = chord_slat*flip(array_slat);
offset_x_slat = gap_x_slat2main + min(nacapane_main.x_pane);                % FOIL 2 X-OFFSET
offset_y_slat = gap_y_slat2main + nacapane_main.y_pane(1);
nacapane_slat.x_pane = array_slat(:, 1) + offset_x_slat;                    % DEFINING THE STRUCT OF POINTS - X
nacapane_slat.y_pane = array_slat(:, 2) + offset_y_slat;                    % DEFINING THE STRUCT OF POINTS - Y

N_slat = length(nacapane_slat.x_pane) - 1;                                  % NUMBER OF PANELS

%% CREATING THE IMAGE OF THE FOILS
nacapane_slat.x_pane_m = array_slat(:, 1) + offset_x_slat;                  % DEFINING THE STRUCT OF POINTS - X (image)
nacapane_slat.y_pane_m = -(array_slat(:, 2) + offset_y_slat);               % DEFINING THE STRUCT OF POINTS - Y (image)
nacapane_main.x_pane_m = array_main(:, 1) ;                                 % DEFINING THE STRUCT OF POINTS - X (image)
nacapane_main.y_pane_m = -(array_main(:, 2) + h_main);                      % DEFINING THE STRUCT OF POINTS - Y (image)
nacapane_flap.x_pane_m = array_flap(:, 1) + offset_x_flap;                  % DEFINING THE STRUCT OF POINTS - X (image)
nacapane_flap.y_pane_m = -(array_flap(:, 2) + offset_y_flap);               % DEFINING THE STRUCT OF POINTS - Y (image)







