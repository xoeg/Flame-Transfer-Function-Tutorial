function [f,FTF] = fn_FTF_exp(path)
    % Function that loops through all the experimental data to compute the
    % FTF. 
    
    % List all of the files in path --------------------------------------
    mylist = ls([path,'\*.mat']);
    % Generate repeating structures --------------------------------------
    WS = load([path,'\',mylist(1,:)]);
    % Generate time vector:
    fs = WS.Params.fs;
    t = 0:(1/fs):((length(WS.P1)-1)/fs);
    % Update Signals structure
    Signals.t = t;
    % Update Measurement structure
    Measurement.Fs = WS.Params.fs;
    Measurement.Mic_Pos = WS.Params.Microphonepos;
    % Geometric and mean structures:
    dp = 0.019; % pipe diameter (mm)
    dr = 0.005; % rod diameter (mm)
    db = 0.013; % bluff body diameter (mm)
    A1 = pi/4 * (dp^2 - dr^2);
    A2 = pi/4 * (dp^2 - db^2);
    % Update Geom structure
    Geom.A1 = A1;
    Geom.A2 = A2;
    % Update the Mean structure
    Mean.c1   = WS.Params.c;
    Mean.rho1 = WS.Params.rho;
    % Here I will assume that the density and the speed of sound do not change
    % with the area contraction:
    Mean.c2   = WS.Params.c;
    Mean.rho2 = WS.Params.rho;
    % But the mean velocity does change:
    u2 = WS.Params.U;
    u1 = u2*A2/A1;
    Mean.u2 = u2; 
    Mean.u1 = u1;
    % Loop through all dataset:
    FTF = zeros(1,size(mylist,1));
    f   = zeros(1,size(mylist,1));
    parfor j = 1:size(mylist,1)
        [f_s,FTF_s] = fun_FTF(j,mylist,path,Signals,Mean,Geom,Measurement);
        FTF(j) = FTF_s;
        f(j) = f_s;
    end
    % Order
    [f,ix] = unique(f);
    FTF = FTF(ix);
    % Plot ---------------------------------------------------------------
    figure
    subplot(2,1,1)
    hold on
    plot(f,abs(FTF),'o')
    xlabel('$f$ (Hz)','Interpreter','Latex')
    ylabel('Gain','Interpreter','Latex')
    ax = gca;
    ax.TickLabelInterpreter = 'Latex';
    box on
    grid on
    ylim([0,max(abs(FTF))])
    
    subplot(2,1,2)
    hold on
    plot(f,unwrap(angle(FTF)),'o')
    xlabel('$f$ (Hz)','Interpreter','Latex')
    ylabel('phase','Interpreter','Latex')
    box on
    grid on
    ax = gca;
    ax.TickLabelInterpreter = 'Latex';
end

function [f,FTF] = fun_FTF(j,mylist,path,Signals,Mean,Geom,Measurement)
    % Function that computes the flame transfer function:   
    WS = load([path,'\',mylist(j,:)]);
    % Generate Signals structure:
    Signals.Forcing = WS.Pref;
    Signals.P_mic = [WS.P1;WS.P2;WS.P3];
    % Multi-Microphone-Method ----------------------------------------
    [F,G,omega] = fn_MMM(Mean,Signals,Measurement);
    % Acoustic field Reconstruction:
    x0 = 0;
    [ph1,uh1] = fn_acoustic_field(Mean,F,G,omega,x0);
    % Compute F and G after the area expansion:
    [F2,G2] = fn_Area_Change(Mean,Geom,ph1,uh1,omega,x0);
    % Compute uh after the area expansion:
    [~,uh] = fn_acoustic_field(Mean,F2,G2,omega,x0);
    % Spectral estimate of PMT signal --------------------------------
    sref = WS.Pref;
    Q    = WS.Q;
    % Remove the mean value of the PMT:
    Q_mean = mean(Q);
    Q = Q - Q_mean;
    % Compute the spectral amplitudes
    t = Signals.t;
    [f,Qh] = fn_spectral_est(t,Q,sref,'welch');
    [Qh_max,ix] = max(Qh);
    freq = f(ix);
    % Compute FTF ---------------------------------------------------
    FTF = (Qh_max/Q_mean)/(uh/Mean.u2);
    f = freq;
end