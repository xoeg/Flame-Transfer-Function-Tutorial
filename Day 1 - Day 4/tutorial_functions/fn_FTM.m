function [F11,F21,FTF] = fn_FTM(Mean,F22)
    % Compute some of the elements of the FTM
    
    % Retrieve mean values:
    rho1 = Mean.rho1;
    rho2 = Mean.rho2;
    c1   = Mean.c1;
    c2   = Mean.c2;
    M1   = Mean.M1;
    gamma = Mean.gamma;
    T2 = Mean.T2;
    T1 = Mean.T1;
    % ====================================================================
    % Compute the following coefficients of the FTM obtained via the
    % Rankine-Hugoniot conditions Eq. (8)
    %   a) F11
    %   b) F21
    %   c) FTF via the element F22


end
