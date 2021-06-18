function [FTF,freqs] = fn_FTF(u1, T1, p2, Qbar, Forcing, Measurement, ...
    Opts)
    % Function that computes the FTF by forcing the flame at different
    % frequencies.     
    
    % Read options:
    fmin = Opts.fmin;
    fmax = Opts.fmax;
    N    = Opts.N;
    % Initialize:
    freqs = linspace(fmin,fmax,N);
    FTF = zeros(1,N);
    % Here you can paralelize the code if you want to make it faster:
    for j = 1:N
        FTF(j) = Flame_Transfer_Function(u1, T1, p2, Qbar, Forcing,...
                 Measurement, freqs, j);
    end
end

function [FTF] = Flame_Transfer_Function(u1, T1, p2, Qbar, Forcing,...
                 Measurement, freqs, j)
    % Computes the Flame transfer function at a single frequency
    % Set the forcing frequency
    Forcing.Frequency = freqs(j);
    %  Run the simulation
    [Mean,~,Signals,~] = Run_Simulation(u1,T1,p2,Qbar,Forcing,Measurement);
    % Read the signals:
    t    = Signals.t;
    sref = Signals.Forcing;
    Q    = Signals.PMT;
    % ====================================================================
    % Cropping the signals -----------------------------------------------
    % Remove the transient part of the PMT signal:
    sref = sref(t > 1);
    Q = Q(t > 1);
    t = t(t > 1);
    % Spectral estimation ------------------------------------------------
    % Remove mean and perform spectral estimation on the rate of heat
    % release:
    % Remove the mean value of the PMT signal:
    Qmean = mean(Q);
    Q = Q - Qmean;
    % Compute the spectral amplitudes:
    [f,Qh] = fn_spectral_est(t,Q,sref,'welch');
    % Find the peak
    [Qh_max,ix] = max(Qh);
    fprintf(['f = ',num2str(f(ix)),' Hz \n'])
    % Multi-Microphone-Method --------------------------------------------
    % Reconstruct the velocity at the reference point:
    [F,G,omega] = fn_MMM(Mean,Signals,Measurement);
    % Signal Reconstruction:
    x_ref = 0;
    [~,uh] = fn_acoustic_field(Mean,F,G,omega,x_ref);
    umean = Mean.u1;
    % Forcing amplitude --------------------------------------------------
    % Compute forcing amplitude and if it is larger than a threshold, say
    % 7%, reduce the forcing level and repeat the simulation.
    A = abs(uh)/umean;
    fprintf(['Forcing Amplitude: ',num2str(A),' \n'])
    if A > 0.07
        V_current = Forcing.Voltage;
        V_new = 0.05*V_current/A;
        Forcing.Voltage = V_new;
        % Recurse:
        [FTF] = Flame_Transfer_Function(u1, T1, p2, Qbar, Forcing,...
                 Measurement, freqs, j);
    else
        % Flame transfer function:
        FTF = (Qh_max/Qmean)/(uh/umean);
    end
end