function [] = fn_plot_acoustic_field(Geom,Mean,Fu,Gu,Fd,Gd,omega,clr)
    % Function that plots the full acoustic field:
    % Retrieve variables:
    Lu = Geom.Lu;
    Lb = Geom.Lb;
    xa = linspace(-Lu,0,1000);
    xb = linspace(0.00001,Lb,1000);
    % Generate acoustic field:
    [pha,uha] = fn_acoustic_field(Mean,Fu,Gu,omega,xa);
    [phb,uhb] = fn_acoustic_field(Mean,Fd,Gd,omega,xb);
    pmax = max(abs([pha,phb]));
    umax = max(abs([uha,uhb]));
    % Plot
    subplot(2,1,1)
    hold on
    plot([xa,xb],abs([pha,phb])/pmax,clr)
    xlabel('$x$ [m]','Interpreter','Latex')
    ylabel('$|\hat{p}|/|\hat{p}_{max}|$','Interpreter','Latex')
    ax = gca;
    ax.TickLabelInterpreter = 'Latex';
    box on
    grid on
    xlim([min(xa) max(xb)])
    
    subplot(2,1,2)
    hold on
    plot([xa,xb],abs([uha,uhb])/umax,clr)
    xlabel('$x$ [m]','Interpreter','Latex')
    ylabel('$|\hat{u}|/|\hat{u}_{max}|$','Interpreter','Latex')
    ax = gca;
    ax.TickLabelInterpreter = 'Latex';
    box on
    grid on
    xlim([min(xa) max(xb)])

end