function [V] = fn_Fit_Distr(omegas,FTF,FB,flag_plot)
    % Fits a a Flame transfer function to two distributed time lag
    % models following AEsoy 2020 CNF "Scaling and prediction..."
    
    %% Setup
    % Fitting options
    optsLS = optimoptions('lsqcurvefit');
    optsLS.Display = 'none';
    % Angular frequencies
    om_end = omegas(end);
    oms = linspace(0,om_end,1000);
    % Separate into gain and phase:
    Gain  = abs(FTF);
    phase = unwrap(angle(FTF));
    %% Step 1: First distribution: =======================================
    % Time delay ---------------------------------------------------------
    % Compute the initial time delay from phase slope: 
    % phase = omega*tau + phi ==> phase = [omegas 1]*[tau;phi]
    [tau_0,phi_0] = LinFit(omegas,-phase);
    % Gain parameters ----------------------------------------------------
    % For the initial values make two sweeps on the parameter space, one
    % for the beta parameter and another for the sigma parameter and look
    % for the solution with the smallest residual norm.
    fun = @(x,w) abs(fn_FTF([exp(0.5*(x(1)*x(2))^2),x,0,0],[0,0,0,0],w));
    N = 10;
    [B,S] = meshgrid(linspace(0,om_end,N),linspace(0,3/om_end,N));
    X0 = [B(:),S(:)];
    x = zeros(size(X0)); Res = zeros(length(X0),1);
    % Plot and loop through initial guesses
    if flag_plot
        Plot_Initial_Guesses_1(X0,omegas,tau_0,phi_0,FTF,N)
    end
    lb = zeros(2,1);
    ub = inf(2,1);
    parfor k = 1:N^2
        [X,R]  = lsqcurvefit(fun,X0(k,:),omegas,Gain,lb,ub,optsLS);
        x(k,:) = X;     Res(k) = R;
    end
    [~,ix] = min(Res); x = x(ix,:);
    % Initial guesses for first distribution:
    v1 = [exp(0.5*(x(1)*x(2))^2),x,tau_0,phi_0];
    if flag_plot
        Fg = figure;
        Plot_Transfer_Function(omegas,FTF,'o',Fg)
        Plot_Transfer_Function(oms,fn_FTF(v1,[0 0 0 0],oms),'--',Fg)
    end
    %% Step 2: Second distribution: =======================================
    % For the second distribution we first analyse the modulations in the
    % gain and the phase. We extract the modulations by removing the first
    % fit, then we try fitting that function to a sine wave and extract the
    % wavelength. We use that wavelength as an initial fit to compute the
    % second time delay:
    FTF_1 = fn_FTF(v1,[0,0,0,0],omegas);
    % Compute the gain difference:
    DG  = Gain - abs(FTF_1);
    % Cut the signal in the frequency band:
    DGr = DG(omegas<=FB(2));
    omb = omegas(omegas<=FB(2));
    DGr = DGr(omb>=FB(1));
    omb = omb(omb>=FB(1));
    % Fit a 1 term fourier series an extract the wave length:
    options = fitoptions('fourier1');
    Lambda = omb(end) - omb(1);
    mx = 2*max(abs(DGr));
    options.Lower = [-1,-mx,-mx,2*pi/Lambda];
    options.Upper = [ 1, mx, mx,6*pi/Lambda];
    ft = fit(omb.',DGr.','fourier1',options);
    Cn = coeffvalues(ft);
    Lambda = 2*pi/Cn(4);
    if flag_plot
       figure; 
       hold on
       plot(omb.',DGr.','o')
       plot(ft)
    end
    % Compute the approximation of the second time delay:
    tau_02 = tau_0 + 2*pi/Lambda;
    % Set up initial guesses. We are going to force the second distribution
    % to be a band pass filter by letting the second distribution to be
    % small at zero frequency.
    tol = 0.01;
    N = 5;
    Betas = linspace(omb(1),omb(end),N);
    Gs = 4*norm(Cn(2:3));
    Sigmas = sqrt(-2*log(tol./Gs)./Betas.^2); 
    [B,S] = deal(zeros(N,N));
    for m = 1:N
        sx1 = log10(Sigmas(m));
        [Bx,Sx] = meshgrid(Betas(m),logspace(sx1,sx1+1,N));
        B(:,m) = Bx;
        S(:,m) = Sx;
    end
    X0 = [Gs*ones(N^2,1),B(:),S(:),tau_02*ones(N^2,1)];
    x = zeros(size(X0)); Res = zeros(length(X0),1);
    % Plot initial guesses
    if flag_plot
        Plot_Initial_Guesses_2(X0,omegas,FTF)
    end
    % Fit the Gain to obtain values of the second distribution:
    fun = @(x,omegas) abs(fn_FTF(v1,x,omegas));
    lb = [0,omb(1),min(min(S)),0];
    ub = [v1(1),inf,max(max(S)),inf];
    parfor k = 1:N^2
        [X,R]  = lsqcurvefit(fun,X0(k,:),omegas,Gain,lb,ub,optsLS);
        x(k,:) = X;
        Res(k) = R;
    end
    [~,ix] = min(Res); v2 = x(ix,:);    
    if flag_plot
        Fg = figure;
        Plot_Transfer_Function(omegas,FTF,'o',Fg)
        Plot_Transfer_Function(oms,fn_FTF(v1,v2,oms),'-',Fg)
        Plot_Transfer_Function(oms,fn_FTF([0,v1(2:5)],v2,oms),'-',Fg)
        % Plot family of results
%         for k = 1:N^2
%             vx = x(k,:);
%             Plot_Transfer_Function(oms,fn_FTF([0,v1(2:5)],vx,oms),'--',Fg)
%             Plot_Transfer_Function(oms,fn_FTF(v1,vx,oms),'--',Fg)
%             kkk = 999;
%         end
    end
    %% Step 3: Full Fit ===================================================
    % Constrain the gain to be 1 at zero frequency such that:
    % g1 = exp(0.5*beta1^2*sigma1^2)*(1 - g2*exp(-0.5*beta2^2*sigma2^2))
    G1 = @(x) exp(0.5*x(1)^2*x(2)^2)*(1-x(5)*exp(-0.5*x(6)^2*x(7)^2));
    fun = @(x,omegas) Split_FTF([G1(x),x(1:4)],x(5:8),omegas);
    FS = [real(FTF).',imag(FTF).'];
    x0 = [v1(2:5),v2];
    lb = zeros(1,8); lb(4) = -pi; lb(6) = omegas(1)/2;
    ub = inf(1,8); ub(4) = pi; ub(5) = 2*v2(1); 
    V = lsqcurvefit(fun,x0,omegas,FS,lb,ub,optsLS);
    G1 = exp(0.5*V(1)^2*V(2)^2)*(1-V(5)*exp(-0.5*V(6)^2*V(7)^2));
    % Outputs ============================================================
    V = [G1 V];
    % Plots ==============================================================
    if flag_plot
       Fg = figure;
       Plot_Transfer_Function(omegas,FTF,'o',Fg)
       Plot_Transfer_Function(oms,fn_FTF(v1,[0,0,0,0],oms),'r--',Fg)
       Plot_Transfer_Function(oms,fn_FTF([0,0,0,0,0],v2,oms),'k--',Fg)
       Plot_Transfer_Function(oms,fn_FTF(V(1:5),[0,0,0,0],oms),'-',Fg)
       Plot_Transfer_Function(oms,fn_FTF([0,0,0,0,V(5)],V(6:9),oms),'-',Fg)
       Plot_Transfer_Function(oms,fn_FTF(V(1:5),V(6:9),oms),'.-',Fg)
    end

    
end
%% Functions
function [m,c] = LinFit(x,y)
    % Function that performs a linear fit: 
    %                      y = m*x + c,
    % where y and x are given and are used to form the matrix:
    %                        A*X = b
    % with A = [x,1] X = [m;c] and b = y. The least squares procedure is
    % performed via QR decomposition.
    
    % Check that x and y are column vectors:
    [s,~] = size(x);
    if s == 1
        x = x.';
    end
    [s,~] = size(y);
    if s == 1
        y = y.';
    end
    % Linear regression
    A = [x ones(size(x))];  
    b = y;
    X = A\b;
    m = X(1);    
    c = X(2);
end
%% Transfer function
function FTF = fn_FTF(x,y,omega)
    % Objective flame transfer function
    G_1     = x(1);
    beta_1  = x(2);
    sigma_1 = x(3);
    Tau_1   = x(4);
    phi     = x(5);
    
    G_2     = y(1);
    beta_2  = y(2);
    sigma_2 = y(3);
    Tau_2   = y(4);
    
    % Transfer function:
    FTF = G_1/2*(exp(-0.5*(omega - beta_1).^2*sigma_1^2) + ...
                 exp(-0.5*(omega + beta_1).^2*sigma_1^2)) ...
                 .* exp(-1i*(omega*Tau_1 + phi)) + ...
          G_2/2*(exp(-0.5*(omega - beta_2).^2*sigma_2^2) + ...
                 exp(-0.5*(omega + beta_2).^2*sigma_2^2)) ...
                 .* exp(-1i*(omega*Tau_2 + phi));
end
function FSP = Split_FTF(x,y,omega)
    % Function that computes a FTF and splits it into real and imaginary
    % parts
    FTF = fn_FTF(x,y,omega);
    FSP(:,1) = real(FTF);
    FSP(:,2) = imag(FTF);
end
%% Plots
function [] = Plot_Transfer_Function(omegas,FTF,Mark,Fig)
    if nargin == 4
        figure(Fig)
    else
        figure
    end
    subplot(2,1,1)
    hold on
    plot(omegas,abs(FTF),Mark)
    box on
    grid on
    xlabel('\omega')
    ylabel('Gain')
    subplot(2,1,2)
    hold on
    plot(omegas,unwrap(angle(FTF)),Mark)
    box on
    grid on
    xlabel('\omega')
    ylabel('phase')
end
function [] = Plot_Initial_Guesses_1(X0,omegas,tau_0,phi_0,FTF,N)
    oms = linspace(0,omegas(end),1000);
    k = 1;
    figure
    for j = 1:N^2
        g1_0 = exp(0.5*(X0(j,1)*X0(j,2))^2);
        FTF_1 = fn_FTF([g1_0,X0(j,:),tau_0,phi_0],[0,0,0,0],oms);
        subplot(3,4,k)
        hold on
        plot(omegas,abs(FTF),'ro-')
        plot(oms,abs(FTF_1),'k-')
        box on
        grid on
        xlabel('\omega')
        ylabel('FTF')
        if mod(j,10) == 0
            k = k + 1;
        end
    end
end
function [] = Plot_Initial_Guesses_2(X0,omegas,FTF)
    oms = linspace(0,omegas(end),1000);
    k = 1;
    figure
    for j = 1:length(X0)
        FTF_1 = fn_FTF([0 0 0 0 0],X0(j,:),oms);
        subplot(2,3,k)
        hold on
        plot(omegas,abs(FTF),'ro-')
        plot(oms,abs(FTF_1),'k-')
        box on
        grid on
        xlabel('\omega')
        ylabel('FTF')
        if mod(j,5) == 0
            k = k + 1;
        end
    end
end