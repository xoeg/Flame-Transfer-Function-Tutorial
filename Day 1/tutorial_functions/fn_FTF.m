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


    % Spectral estimation ------------------------------------------------
    % Remove mean and perform spectral estimation on the rate of heat
    % release:


    % Multi-Microphone-Method --------------------------------------------
    % Reconstruct the velocity at the reference point:


    % Forcing amplitude --------------------------------------------------
    % Compute forcing amplitude and if it is larger than a threshold, say
    % 7%, reduce the forcing level and repeat the simulation.


end
