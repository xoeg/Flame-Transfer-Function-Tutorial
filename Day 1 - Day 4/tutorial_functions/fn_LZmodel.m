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

end
