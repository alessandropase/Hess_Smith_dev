function u_tot = TotalVelocity_groundX3 (q_1, q_2, q_3, gamma_1, ...
                                         gamma_2, gamma_3, u_inf, x, y, ...
                                         nacapane_1, nacapane_2, ...
                                         nacapane_3, N1, N2, N3)

u_tot_1 = 0;
u_tot_2 = 0;
u_tot_3 = 0;

% 1

for i = 1 : N1
    beta = atan2((nacapane_1.y_pane(i+1) - nacapane_1.y_pane(i)) , ...
                 (nacapane_1.x_pane(i+1) - nacapane_1.x_pane(i)));
    Q = [ cos(beta), sin(beta);
         -sin(beta), cos(beta)
        ];

    Estremo_1 = [nacapane_1.x_pane(i); nacapane_1.y_pane(i)];
    Estremo_2 = [nacapane_1.x_pane(i+1); nacapane_1.y_pane(i+1)];
    
    
    u_sorg = ViSorgente([x;y], Estremo_1, Estremo_2, Q', Q);
    u_vort = ViVortice([x;y], Estremo_1, Estremo_2, Q', Q);

    
    u_tot_1 = u_tot_1 + u_sorg * q_1(i) + gamma_1 * u_vort;
end

% 2

for i = 1 : N2
    beta = atan2((nacapane_2.y_pane(i+1) - nacapane_2.y_pane(i)) , ...
                 (nacapane_2.x_pane(i+1) - nacapane_2.x_pane(i)));
    Q = [ cos(beta), sin(beta);
         -sin(beta), cos(beta)
        ];

    Estremo_1 = [nacapane_2.x_pane(i); nacapane_2.y_pane(i)];
    Estremo_2 = [nacapane_2.x_pane(i+1); nacapane_2.y_pane(i+1)];
    
    
    u_sorg = ViSorgente([x;y], Estremo_1, Estremo_2, Q', Q);
    u_vort = ViVortice([x;y], Estremo_1, Estremo_2, Q', Q);
    
    
    u_tot_2 = u_tot_2 + u_sorg * q_2(i) + gamma_2 * u_vort;
end

% 3
for i = 1 : N3
    beta = atan2((nacapane_3.y_pane(i+1) - nacapane_3.y_pane(i)) , ...
                 (nacapane_3.x_pane(i+1) - nacapane_3.x_pane(i)));
    Q = [ cos(beta), sin(beta);
         -sin(beta), cos(beta)
        ];

    Estremo_1 = [nacapane_3.x_pane(i); nacapane_3.y_pane(i)];
    Estremo_2 = [nacapane_3.x_pane(i+1); nacapane_3.y_pane(i+1)];
    
    
    u_sorg = ViSorgente([x;y], Estremo_1, Estremo_2, Q', Q);
    u_vort = ViVortice([x;y], Estremo_1, Estremo_2, Q', Q);
    
    
    u_tot_3 = u_tot_3 + u_sorg * q_3(i) + gamma_3 * u_vort;
end

% tot
u_tot = u_tot_1 + u_tot_2 + u_tot_3 + u_inf;
