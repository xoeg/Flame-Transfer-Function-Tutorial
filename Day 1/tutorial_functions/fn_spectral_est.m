function [f,X,PSD] = fn_spectral_est(t,x,xref,type,Opts)
    % Function that computes the spectral estimation of a real signal,
    % cross correlated with a reference signal, using a direct windowed
    % method or welch method. 
    % Inputs:
    %           t - time vector
    %           x - real signal (column vector)
    %           xref - reference signal
    %           type - Can be two options
    %                   'direct' - simple fft
    %                   'welch'  - welch method
    %           Opts - Optional structure with fields:
    %                   N - Zero padding length
    %                   M - Number of segments for welch mehtod
    %                   p - Overlapping percentage 
    % Outputs:
    %           f = vector of frequencies
    %           X = Fourier transformed signal (zero padded, windowed and
    %           or averaged) 
    %           PSD = Power spectral density

    % Sanity checks ------------------------------------------------------
    % Ensure x is a real signal:
    if sum(imag(x)) ~= 0
        error('signal x is not real')
    end
    % Ensure that x is a single column vector
    sz = size(x);
    if ~any(sz == 1)
        error('Signal of unexpected size')
    elseif sz(2) > 1
        x = x.';
    end
    % Ensure that xref is a single column vector
    if ~isempty(xref)
        sz = size(xref);
        if ~any(sz == 1)
            error('Signal of unexpected size')
        elseif sz(2) > 1
            xref = xref.';
        end
    end

    
    % Basic definitions --------------------------------------------------
    L  = length(t);           % Length of signal
    Ts = t(2) - t(1);         % Sampling period
    Fs = 1/Ts;                % Sampling frequency
    fmin = 0.01;              % Minimum resolved frequency (Hz)
    N = 2^nextpow2(Fs/fmin);  % Zero padded length
    % Options ------------------------------------------------------------
    % For welch method, define the default number of segments "M" and the
    % percentage of overlapping "p":
    M = 8;
    p = 0.5;
    % Check for optional inputs:
    if nargin == 5
        % Change default zero padding length:
        if isfield(Opts,'N')
            N = Opts.N;
        end
        if isfield(Opts,'M')
            M = Opts.M;
        end
        if isfield(Opts,'p')
            p = Opts.p;
        end
    end
    % Since the signal is real we are only interested in the single sided
    % spectra:
    % Vector of frequencies:
    f = Fs*(0:(N/2))/N;
    % Compute amplitude and PSD
    switch type
        case 'direct'
            % Windowing --------------------------------------------------
            % Hanning window
            w = hann(L);
            ACF = 1/mean(w);    % Amplitude correction factor
            ECF = 1/rms(w);     % Energy correction factor
            % FFT --------------------------------------------------------
            X = fft(w.*x,N);
            % Cross correlation
            if ~isempty(xref)
                Xref = fft(w.*xref,N);
                XXref    = X    .* conj(Xref);
                XrefXref = Xref .* conj(Xref);
                X = XXref ./ sqrt(XrefXref);
            end
            % PSD estimation (single sided) ------------------------------
            PSD = abs(X * ECF).^2 * Ts/L;
            PSD = PSD(1:N/2+1);
            PSD(2:end-1) = 2*PSD(2:end-1);
            % Amplitude estimation (single sided) ------------------------
            X = X * ACF/L;
            X = X(1:N/2+1);
            X(2:end-1) = 2*X(2:end-1);            
        case 'welch'
            % Splitting the signal into "M" segments with "p" overlap. 
            Lm = ceil(L/(M*(1-p)));     % Length of each segment 
            Lo = floor(Lm*p);           % Overlap length
            % Splitting signal:
            xs = buffer(x,Lm,Lo,'nodelay');
            if ~isempty(xref)
                xrefs = buffer(xref,Lm,Lo,'nodelay');
            end
            % Windowing --------------------------------------------------
            % Hanning window
            w = hann(Lm);
            ACF = 1/mean(w);    % Amplitude correction factor
            ECF = 1/rms(w);     % Energy correction factor
            % Loop each segment of the signal:
            M = size(xs,2);
            PSD = zeros(N/2+1,1);
            X = zeros(N/2+1,1);
            for j = 1:M
                % FFT ----------------------------------------------------
                Xs = fft(w.*xs(:,j),N);
                % Cross correlation
                if ~isempty(xref)
                    Xref = fft(w.*xrefs(:,j),N);
                else
                    % If there is no reference signal, auto correlate to
                    % preserve the amplitude at the expense of losing the
                    % phase of the signal.
                    Xref = fft(w.*xs(:,j),N);
                end
                XXref    = Xs   .* conj(Xref);
                XrefXref = Xref .* conj(Xref);
                Xs = XXref ./ sqrt(XrefXref);
                % PSD estimation (single sided) --------------------------
                Ps = abs(Xs * ECF).^2 * Ts/Lm;
                Ps = Ps(1:N/2+1);
                Ps(2:end-1) = 2*Ps(2:end-1);
                PSD = PSD + Ps;
                % Amplitude estimation (single sided) --------------------
                Xs = Xs * ACF/Lm;
                Xs = Xs(1:N/2+1);
                Xs(2:end-1) = 2*Xs(2:end-1);
                X = X + Xs;
            end
            % Average:
            PSD = PSD/M;
            X = X/M;
    end
end