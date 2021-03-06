function [F,G,omega] = fn_MMM(Mean,Signals,Measurement)
    % Function that computes the values of the wave amplitudes from the
    % Signals
    
    % Retrieve signals and postitions:
    t    = Signals.t;
    if isfield(Signals,'Forcing')
        xref = Signals.Forcing;
    elseif isfield(Signals,'Forcing_Downstream')
        xref = Signals.Forcing_Downstream;
    end
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
    indx1 = find(t == 1);
    t     = t(indx1:end);
    xref  = xref(indx1:end);
    pp    = pp(:,indx1:end);
    % Spectral estimation ------------------------------------------------
    % Using the function fn_spectral_est compute the values of p-hat (ph)
    % and omega for each pressure time series. The pressure time series are
    % given by the matrix "pp". The jth column of "pp" i.e., pp(j,:)
    % corresponds to the pressure time series at location x_j = x_mic(j).
    [ph,omega] = deal(zeros(size(pp,1),1));
    for j = 1:size(pp,1)
        [f,Ph] = fn_spectral_est(t,pp(j,:),xref,'welch');
        % Find peak:
        [ph(j),ix] = max(Ph);
        omega(j) = 2*pi*f(ix);
    end
    % Take the average of the frequencies
    omega = mean(omega);
    % Find the wave amplitudes -------------------------------------------
    % Using Eq. (4) find the wave amplotudes F and G
    % Compute wave numbers:
    kp = - omega/(c + u);
    km =   omega/(c - u);
    % Find F and G (Riemman invariants) via least squares. 
    A = [exp(1i*kp*x_mic) exp(1i*km*x_mic)];
    b = ph;
    % Least squares via qr decomposition:
    x = A\b;
    F = x(1);
    G = x(2);
end