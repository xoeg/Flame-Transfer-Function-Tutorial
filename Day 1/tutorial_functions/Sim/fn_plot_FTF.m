function [] = fn_plot_FTF(f,FTF,Mean,Geom)
    % Function to plot the Flame transfer function


    % Analytical Flame Transfer Function ---------------------------------
    % Read parameters:
    b   = Geom.b;
    a   = Geom.a;
    Su  = Mean.Su;
    ug  = Mean.u1;
    Tau = Mean.Tau;
    % Compute analytic FTF
    fa = linspace(0.001,400,1000);
    omega = 2 * pi * fa;
    Omega = omega * (b-a)/(Su * sqrt(1 - Su^2/ug^2));
    FTFa = 2./(1i*Omega*(a + b)) .* (a - b*exp(-1i*Omega) + ...
           (b-a)./(1i*Omega).*(1-exp(-1i*Omega)));
    FTF_analytic = FTFa .* exp(-1i*omega*Tau);
    % Plot FTFs and compare ----------------------------------------------
    figure
    subplot(2,1,1)
    hold on
    plot(f,abs(FTF),'o')
    plot(fa,abs(FTF_analytic))
    xlabel('$f$ (Hz)','Interpreter','Latex')
    ylabel('Gain','Interpreter','Latex')
    legend('Measured','Analytical','Location','northeast')
    ax = gca;
    ax.TickLabelInterpreter = 'Latex';
    box on
    grid on
    
    % Unwrap the experimental phase according to the analytical:
    phase_analytic = unwrap(angle(FTF_analytic));
    phase_exp = unwrap(angle(FTF));
    for j = 1:length(f)
        [~,ix] = min(abs(fa - f(j)));
        pha = phase_analytic(ix);
        phe = phase_exp(j);
        % Difference
        dp = (pha - phe)/(2*pi);
        % Update:
        phase_exp(j) = phe + sign(dp)*floor(abs(dp))*2*pi;
    end
    
    subplot(2,1,2)
    hold on
    plot(f,phase_exp,'o')
    plot(fa,phase_analytic)
    xlabel('$f$ (Hz)','Interpreter','Latex')
    ylabel('phase','Interpreter','Latex')
    box on
    grid on
    ax = gca;
    ax.TickLabelInterpreter = 'Latex';

end