function u_tot = TotalVelocity_ground (q, gamma, u_inf, x, y, x_pane, y_pane, x_pane_m, y_pane_m)

N = length(x_pane)-1;
u_tot = u_inf;

for i = 1 : N
    beta = atan2((y_pane(i+1) - y_pane(i)) , ...
                 (x_pane(i+1) - x_pane(i)));
    Q = [ cos(beta), sin(beta);
         -sin(beta), cos(beta)
        ];

    Estremo_1 = [x_pane(i); y_pane(i)];
    Estremo_2 = [x_pane(i+1); y_pane(i+1)];
    
    Estremo_1_m = [x_pane_m(i); y_pane_m(i)];               
    Estremo_2_m = [x_pane_m(i+1); y_pane_m(i+1)];           
    
    beta_m = atan2((Estremo_2_m(2) - Estremo_1_m(2)) , ...
                   (Estremo_2_m(1) - Estremo_1_m(1)));    
    Q_m = [ cos(beta_m), sin(beta_m);
           -sin(beta_m), cos(beta_m)
          ]; 

    
    u_sorg = ViSorgente([x;y], Estremo_1, Estremo_2, Q', Q);
    u_vort = ViVortice([x;y], Estremo_1, Estremo_2, Q', Q);
    u_sorg_m = ViSorgente([x;y], Estremo_1_m, Estremo_2_m, Q_m', Q_m);
    u_vort_m = ViVortice([x;y], Estremo_1_m, Estremo_2_m, Q_m', Q_m);


    
    u_tot = u_tot + (u_sorg+u_sorg_m)*q(i) + gamma*(u_vort-u_vort_m);
end
