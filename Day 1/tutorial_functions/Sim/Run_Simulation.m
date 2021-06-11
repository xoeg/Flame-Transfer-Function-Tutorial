function [Mean,Geom,Signals,Sim] = Run_Simulation(u1, T1, p2, Qbar, ...
     Forcing, Measurement, order)
    % Wrapping function to run a simulation

    % Set up =============================================================
    % Set up measurements:
    % Hotwire:
    if strcmp(Measurement.HWA,'ON') && length(Measurement.HWA_Pos) > 1
        error('Our Lab has only one hot wire anemometer.')
    end
    % Define spatial discretization order:
    if nargin == 7 && order == 1
        order = 1;
    else
        order = 2;
    end
    % Load the geometry
    Geom = fn_define_geometry(order);
    % Compute the mean flow:
    Mean = fn_mean_flow(u1,T1,p2,Qbar,Geom);  
    % Set up Forcing:
    if strcmp(Forcing.Speaker,'ON')
        Geom.Speaker = true;
        % Set up Forcing:
        Mean.Force.V = Forcing.Voltage * 1e-2;
        Mean.Force.omega = Forcing.Frequency * 2*pi;
    else
        Geom.Speaker = false;
        % Set up Forcing:
        Mean.Force.V = 0;
        Mean.Force.omega = 0;
    end
    
    % Time marching the simulation =======================================
    dt    = 1/10007;               % Time step of simulation
    t_end = Measurement.Time;      % Length of simulation
    f_sim = 1/dt;                  % Simulation sampling frequency
    if Measurement.Fs > f_sim
        error(['Lower sampling frequency below: ',num2str(f_sim),' Hz'])
    end
    [ts,Gs,Hs,Xis,Qs] = fn_time_marching(t_end,dt,Geom,Mean);
    % Save simulation values
    Sim.Geom = Geom;
    Sim.Mean = Mean;
    Sim.ts = ts;
    Sim.Qs = Qs;
    Sim.Gs = Gs;
    Sim.Hs = Hs;
    Sim.Xis = Xis;
    Sim.Forcing = Forcing;
    Sim.Measurement = Measurement;
    % Measurements =======================================================
    % All resampling is done assuming linear interpolation
    % Resampling:
    fs = Measurement.Fs;
    t = 0:1/fs:ts(end);
    Signals.t = t;
    % Forcing signal -----------------------------------------------------
    if strcmp(Forcing.Speaker,'ON')
        V     = Mean.Force.V / 1e-2;
        omega = Mean.Force.omega;
        F     = V*cos(omega*t);
        Signals.Forcing = F;
    end
    % Photomultiplier ----------------------------------------------------
    if strcmp(Measurement.PMT,'ON')
        Qrs = interp1(ts,Qs,t);
        Signals.PMT = Qrs;
    end
    % Microphones: -------------------------------------------------------
    if strcmp(Measurement.Mic,'ON')
        x_mic = Measurement.Mic_Pos;
        Acoustic_only = true;
        Lu = Geom.Lu;
        Lb = Geom.Lb;
        P_mic = zeros(length(x_mic),length(t));
        for k = 1:length(x_mic)
            if x_mic(k) < -Lu || x_mic(k) > Lb
                % If location is out of the combustor, send white noise:
                P_mic(k,:) = randn(1,length(t))*10^-3;
            elseif x_mic(k) > 0 && Mean.T2 > 373
                % If location is in the hot section and temperature is
                % above 100 C, fry the microphone after a few sampling
                % periods and send noise. 
                t_burn = (200 + round(rand()*100))/fs;
                for j = 1:length(t)
                    if t(j) < t_burn
                        [~,pp] = fn_modehshapes(t(j),x_mic(k),ts,Gs,Hs,...
                                              Mean,Geom,Acoustic_only);
                        P_mic(k,j) = pp;
                    else
                        P_mic(k,j) = randn();
                    end
                end
            else
                for j = 1:length(t)
                    [~,pp] = fn_modehshapes(t(j),x_mic(k),ts,Gs,Hs,Mean,...
                                            Geom,Acoustic_only);
                    P_mic(k,j) = pp;
                end
            end
        end
        Signals.P_mic = P_mic;
    end
    % Hot wire: ----------------------------------------------------------
    if strcmp(Measurement.HWA,'ON')
        x_hw = Measurement.HWA_Pos;
        Acoustic_only = false;
        Lu = Geom.Lu;
        Lb = Geom.Lb;
        U_hw = zeros(length(x_hw),length(t));
        for k = 1:length(x_hw)
            if x_hw(k) < -Lu || x_hw(k) > Lb
                % If location is out of the combustor, send white noise:
                U_hw(k,:) = randn(1,length(t))*10^-3;
            elseif x_hw(k) > 0 && Mean.T2 > 373
                % If location is in the hot section and temperature is
                % above 100 C, fry the microphone after a few sampling
                % periods and send noise. 
                t_burn = (200 + round(rand()*100))/fs;
                for j = 1:length(t)
                    if t(j) < t_burn
                        [~,pp] = fn_modehshapes(t(j),x_hw(k),ts,Gs,Hs,...
                                              Mean,Geom,Acoustic_only);
                        U_hw(k,j) = pp;
                    else
                        U_hw(k,:) = randn();
                    end
                end
            else
                for j = 1:length(t)
                    [up] = fn_modehshapes(t(j),x_hw(k),ts,Gs,Hs,Mean,...
                                            Geom,Acoustic_only);
                    U_hw(k,j) = up;
                end
            end
        end
        Signals.U_hw = U_hw;
    end
    fprintf('\n Finished \n')
end