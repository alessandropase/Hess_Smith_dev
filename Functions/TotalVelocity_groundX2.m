function u_tot = TotalVelocity_groundX2 (q_1, q_2, gamma_1, gamma_2, u_inf, x, y, nacapane_1, nacapane_2, N1, N2)

u_tot_1 = 0;
u_tot_2 = 0;

for i = 1 : N1
    beta = atan2((nacapane_1.y_pane(i+1) - nacapane_1.y_pane(i)) , ...
                 (nacapane_1.x_pane(i+1) - nacapane_1.x_pane(i)));
    Q = [ cos(beta), sin(beta);
         -sin(beta), cos(beta)
        ];

    Estremo_1 = [nacapane_1.x_pane(i); nacapane_1.y_pane(i)];
    Estremo_2 = [nacapane_1.x_pane(i+1); nacapane_1.y_pane(i+1)];
    
    Estremo_1_m = [nacapane_1.x_pane_m(i); nacapane_1.y_pane_m(i)];               
    Estremo_2_m = [nacapane_1.x_pane_m(i+1); nacapane_1.y_pane_m(i+1)];           
    
    beta_m = atan2((Estremo_2_m(2) - Estremo_1_m(2)) , ...
                   (Estremo_2_m(1) - Estremo_1_m(1)));    
    Q_m = [ cos(beta_m), sin(beta_m);
           -sin(beta_m), cos(beta_m)
          ]; 

    
    u_sorg = ViSorgente([x;y], Estremo_1, Estremo_2, Q', Q);
    u_vort = ViVortice([x;y], Estremo_1, Estremo_2, Q', Q);
    u_sorg_m = ViSorgente([x;y], Estremo_1_m, Estremo_2_m, Q_m', Q_m);
    u_vort_m = ViVortice([x;y], Estremo_1_m, Estremo_2_m, Q_m', Q_m);


    
    u_tot_1 = u_tot_1 + (u_sorg+u_sorg_m)*q_1(i) + gamma_1*(u_vort-u_vort_m);
end

for i = 1 : N2
    beta = atan2((nacapane_2.y_pane(i+1) - nacapane_2.y_pane(i)) , ...
                 (nacapane_2.x_pane(i+1) - nacapane_2.x_pane(i)));
    Q = [ cos(beta), sin(beta);
         -sin(beta), cos(beta)
        ];

    Estremo_1 = [nacapane_2.x_pane(i); nacapane_2.y_pane(i)];
    Estremo_2 = [nacapane_2.x_pane(i+1); nacapane_2.y_pane(i+1)];
    
    Estremo_1_m = [nacapane_2.x_pane_m(i); nacapane_2.y_pane_m(i)];               
    Estremo_2_m = [nacapane_2.x_pane_m(i+1); nacapane_2.y_pane_m(i+1)];           
    
    beta_m = atan2((Estremo_2_m(2) - Estremo_1_m(2)) , ...
                   (Estremo_2_m(1) - Estremo_1_m(1)));    
    Q_m = [ cos(beta_m), sin(beta_m);
           -sin(beta_m), cos(beta_m)
          ]; 

    
    u_sorg = ViSorgente([x;y], Estremo_1, Estremo_2, Q', Q);
    u_vort = ViVortice([x;y], Estremo_1, Estremo_2, Q', Q);
    u_sorg_m = ViSorgente([x;y], Estremo_1_m, Estremo_2_m, Q_m', Q_m);
    u_vort_m = ViVortice([x;y], Estremo_1_m, Estremo_2_m, Q_m', Q_m);


    
    u_tot_2 = u_tot_2 + (u_sorg+u_sorg_m)*q_2(i) + gamma_2*(u_vort-u_vort_m);
end

u_tot = u_tot_1 + u_tot_2 + u_inf;