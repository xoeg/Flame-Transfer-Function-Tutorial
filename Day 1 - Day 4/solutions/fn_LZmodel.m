function [B11,B22,Leq,Zeta] = fn_LZmodel(Mean,Geom)
    % Computes components of the BTM for the L-Z model.
    
    % Retrieve mean values:
    rho1 = Mean.rho1;
    rho2 = Mean.rho2;
    c1 = Mean.c1;
    c2 = Mean.c2;
    M1 = Mean.M1;
    A1 = Geom.A1;
    A2 = Geom.A2;
    % Load values for the B12 element:
    WS = load('BTM_1.mat');
    B12 = WS.B12;
    freqs = WS.freqs;
    omegas = 2*pi*freqs.';
    
    % ===================================================================
    % Compute the following coefficients of the FTM obtained via the
    % L-Z model Eq. (5):
    %   a) B11
    %   b) B22
    %   c) Leq and Zeta. For this, the following link might be helpful:
    % https://uk.mathworks.com/help/optim/ug/fit-model-to-complex-data.html
    % Look at the section: Alternative: Split Real and Imaginary Parts.
    % Assume that Leq and Zeta are real numbers!
    B11 = rho1*c1/(rho2*c2);
    B22 = A1/A2;
    % Define the loss coefficient function:
    fun_LZ = @(L,Z,omegas) rho1*c1/(rho2*c2)*(-1i*omegas*L/c1 - M1*Z);
    % Set the least squares routine
    x0 = [1,1]; % Arbitrary initial values
    obj = @(v,omegas) [real(fun_LZ(v(1),v(2),omegas)) imag(fun_LZ(v(1),v(2),omegas))];
    opts = optimoptions(@lsqcurvefit,'Display','none');
    [v_est] = lsqcurvefit(obj,x0,omegas,[real(B12) imag(B12)],[],[],opts);
    % L-Z estimates:
    Leq  = v_est(1);
    Zeta = v_est(2);
end