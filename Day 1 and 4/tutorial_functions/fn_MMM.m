function [F,G,omega] = fn_MMM(Mean,Signals,Measurement)
    % Function that computes the values of the wave amplitudes from the
    % Signals
    
    % Retrieve signals and postitions:
    t    = Signals.t;
    xref = Signals.Forcing;
    pp   = Signals.P_mic;
    x_mic = Measurement.Mic_Pos;
    % Retrieve mean flow values 
    if all(x_mic <= 0)
        u = Mean.u1;
        c = Mean.c1;
    elseif all(x_mic > 0)
        u = Mean.u2;
        c = Mean.c2;
    else
        error('Check microphone positions')
    end
    % Sanity checks ------------------------------------------------------
    % 1) Ensure x_mic is a column vector:
    if size(x_mic,2) > 1 && size(x_mic,1) == 1
        x_mic = x_mic.';
    else
        error('Expecting a vector')
    end
    % ====================================================================
    % Cropping the signals -----------------------------------------------
    % Remove the transient part of the signals as done in the example:

    
    
    % Spectral estimation ------------------------------------------------
    % Using the function fn_spectral_est compute the values of p-hat (ph)
    % and omega for each pressure time series. The pressure time series are
    % given by the matrix "pp". The jth column of "pp" i.e., pp(j,:)
    % corresponds to the pressure time series at location x_j = x_mic(j).

    
    
    
    
    % Find the wave amplitudes -------------------------------------------
    % Using Eq. (4) find the wave amplotudes F and G

    
    
    
end