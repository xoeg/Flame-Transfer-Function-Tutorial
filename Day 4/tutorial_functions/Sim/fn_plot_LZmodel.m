function [] = fn_plot_LZmodel(Mean,Leq,Zeta)
    
    % Retreive Mean values
    rho2 = Mean.rho2;
    rho1 = Mean.rho1;
    c2   = Mean.c2;
    c1   = Mean.c1;
    M1   = Mean.M1;

    % Load values for the B12 element:
    WS = load('BTM_1.mat');
    B12 = WS.B12;
    freqs = WS.freqs;
    omegas = 2*pi*freqs.';
    % LZ model
    LZ_model = rho1*c1/(rho2*c2)*(-1i*omegas*Leq/c1 - M1*Zeta);
        
    % Plot and compare:
    figure
    subplot(2,1,1)
    hold on
    plot(freqs,abs(B12))
    plot(freqs,abs(LZ_model))
    xlabel('$f$ [Hz]','Interpreter','Latex')
    ylabel('$|B_{12}|$','Interpreter','Latex')
    box on
    grid on
    ax = gca;
    ax.TickLabelInterpreter = 'Latex';
    
    subplot(2,1,2)
    hold on
    plot(freqs,angle(B12))
    plot(freqs,angle(LZ_model))
    xlabel('$f$ [Hz]','Interpreter','Latex')
    ylabel('$\angle(B_{12})$','Interpreter','Latex')
    box on
    grid on
    ax = gca;
    ax.TickLabelInterpreter = 'Latex';

end