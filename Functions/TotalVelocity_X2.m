function u_tot = TotalVelocity_X2 (q_1, q_2, gamma_1, gamma_2, u_inf, x, y, ...
                                  x_pane_1, y_pane_1, x_pane_2, y_pane_2)

N = length(x_pane_1)-1;
u_tot = u_inf;

% u foil 1
for i = 1 : N
    beta = atan2((y_pane_1(i+1) - y_pane_1(i)) , ...
                 (x_pane_1(i+1) - x_pane_1(i)));
    Q = [ cos(beta), sin(beta);
         -sin(beta), cos(beta)
        ];

    Estremo_1 = [x_pane_1(i); y_pane_1(i)];
    Estremo_2 = [x_pane_1(i+1); y_pane_1(i+1)];

    u_sorg = ViSorgente([x;y], Estremo_1, Estremo_2, Q', Q);
    u_vort = ViVortice([x;y], Estremo_1, Estremo_2, Q', Q);
    
    u_tot = u_tot + u_sorg*q_1(i) + gamma_1*u_vort;
end

% u foil 2
for i = 1 : N
    beta = atan2((y_pane_2(i+1) - y_pane_2(i)) , ...
                 (x_pane_2(i+1) - x_pane_2(i)));
    Q = [ cos(beta), sin(beta);
         -sin(beta), cos(beta)
        ];

    Estremo_1 = [x_pane_2(i); y_pane_2(i)];
    Estremo_2 = [x_pane_2(i+1); y_pane_2(i+1)];

    u_sorg = ViSorgente([x;y], Estremo_1, Estremo_2, Q', Q);
    u_vort = ViVortice([x;y], Estremo_1, Estremo_2, Q', Q);
    
    u_tot = u_tot + u_sorg*q_2(i) + gamma_2*u_vort;
end
