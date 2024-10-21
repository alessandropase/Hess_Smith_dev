function Velocity_Plot(Nfx, Nfy, x_lim, y_lim, q_1, q_2, q_3, gamma_1, ...
                       gamma_2, gamma_3, nacapane_1, nacapane_2, ...
                       nacapane_3, N1, N2, N3, u_inf, Field)

x = linspace(x_lim(1),x_lim(2),Nfx);
y = linspace(y_lim(1),y_lim(2),Nfy);
Ufield = zeros(Nfx,Nfy); 
Vfield = zeros(Nfx,Nfy); 
CPfield = zeros(Nfx,Nfy); 

for i = 1:Nfx
    for j = 1:Nfy
        u_aux = TotalVelocity_groundX3 (q_1, q_2, q_3, gamma_1, gamma_2, gamma_3, u_inf, x(i), y(j), nacapane_1, nacapane_2, nacapane_3, N1, N2, N3);
        Vfield(i,j) = u_aux'*[0;1];
        Ufield(i,j) = u_aux'*[1;0];
        CPfield(i,j) = 1 - (norm(u_aux)^2)/(norm(u_inf)^2);
    end
end

[X,Y] = meshgrid(x,y);

if Field == 'U_Field' | Field == 'A_Field'
    figure;
    hold on
    box on
    grid on
    contourf(X',Y', Ufield,100,'LineStyle','None');
    colormap(flipud(hot));
    colorbar;
    c=colorbar('northoutside');
    c.TickLabelInterpreter='latex';
    %colorbar('off')
    plot(nacapane_1.x_pane,nacapane_1.y_pane,'Linewidth',2,'Color','k')
    plot(nacapane_2.x_pane,nacapane_2.y_pane,'Linewidth',2,'Color','k')
    plot(nacapane_3.x_pane,nacapane_3.y_pane,'Linewidth',2,'Color','k')
    view(2);
    tx=xlabel('$x$');
    ty=ylabel('$y$');
    tx.Interpreter='latex';
    ty.Interpreter='latex';
    set(gca,'TickLabelInterpreter', 'latex');
    x0=10;
    y0=10;
    width=1250;
    height=580;
    set(gcf,'position',[x0,y0,width,height])
    ax = gca;
    axis equal
    axis off
    ax.FontSize = 40;
    fill(nacapane_1.x_pane,nacapane_1.y_pane, [0.7, 0.7, 0.7])
    fill(nacapane_2.x_pane,nacapane_2.y_pane, [0.7, 0.7, 0.7])
    fill(nacapane_3.x_pane,nacapane_3.y_pane, [0.7, 0.7, 0.7])
end

if Field == 'V_Field' | Field == 'A_Field'
    figure;
    hold on
    box on
    grid on
    contourf(X',Y',Vfield,100,'LineStyle','None');
    colormap(bluewhitered_aero(256)); 
    colorbar;
    c=colorbar('northoutside');
    c.TickLabelInterpreter='latex';
    %colorbar('off')
    plot(nacapane_1.x_pane,nacapane_1.y_pane,'Linewidth',2,'Color','k')
    plot(nacapane_2.x_pane,nacapane_2.y_pane,'Linewidth',2,'Color','k')
    plot(nacapane_3.x_pane,nacapane_3.y_pane,'Linewidth',2,'Color','k')
    view(2);
    tx=xlabel('$x$');
    ty=ylabel('$y$');
    tx.Interpreter='latex';
    ty.Interpreter='latex';
    set(gca,'TickLabelInterpreter', 'latex');  
    x0=10;
    y0=10;
    width=1250; 
    height=580;
    set(gcf,'position',[x0,y0,width,height])
    ax = gca;
    axis equal
    axis off
    ax.FontSize = 40;
    fill(nacapane_1.x_pane,nacapane_1.y_pane, [0.7, 0.7, 0.7])
    fill(nacapane_2.x_pane,nacapane_2.y_pane, [0.7, 0.7, 0.7])
    fill(nacapane_3.x_pane,nacapane_3.y_pane, [0.7, 0.7, 0.7])
end

if  Field == 'P_Field' | Field == 'A_Field'
    figure;
    hold on
    box on
    grid on
    contourf(X',Y', CPfield,100,'LineStyle','None');
    colormap(flipud(hot));
    clim([-1.5,0.87])
    colorbar;
    c=colorbar('northoutside');
    c.TickLabelInterpreter='latex';
    %colorbar('off')
    plot(nacapane_1.x_pane,nacapane_1.y_pane,'Linewidth',2,'Color','k')
    plot(nacapane_2.x_pane,nacapane_2.y_pane,'Linewidth',2,'Color','k')
    plot(nacapane_3.x_pane,nacapane_3.y_pane,'Linewidth',2,'Color','k')
    view(2);
    tx=xlabel('$x$');
    ty=ylabel('$y$');
    tx.Interpreter='latex';
    ty.Interpreter='latex';
    set(gca,'TickLabelInterpreter', 'latex');  
    x0=10;
    y0=10;
    width=1250; 
    height=580;
    set(gcf,'position',[x0,y0,width,height])
    ax = gca;
    axis equal
    axis off
    ax.FontSize = 40;
    fill(nacapane_1.x_pane,nacapane_1.y_pane, [0.7, 0.7, 0.7])
    fill(nacapane_2.x_pane,nacapane_2.y_pane, [0.7, 0.7, 0.7])
    fill(nacapane_3.x_pane,nacapane_3.y_pane, [0.7, 0.7, 0.7])
end