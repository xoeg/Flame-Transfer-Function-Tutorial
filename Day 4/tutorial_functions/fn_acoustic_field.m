function [ph,uh] = fn_acoustic_field(Mean,F,G,omega,x)
    % Function that recontructs the acoustic field using the wave
    % amplitudes.
    % The function receives as an input the "Mean" flow structure, the wave
    % amplitudes "F", "G", the angular frequency "omega", and a vector of
    % positions "x". It outputs a vector of complex pressure "ph" and a
    % vector of complex velocities "uh".
    
    % Retrieve mean flow values 
    if all(x <= 0)
        rho = Mean.rho1;
        c = Mean.c1;
        u = Mean.u1;
    elseif all(x > 0)
        rho = Mean.rho2;
        c = Mean.c2;
        u = Mean.u2;
    else
        error('Flow properties change in the selected x range.')
    end
    % ====================================================================
    % Reconstruct the acoustic field by writing the equations for p-hat
    % (ph) and u-hat (uh) Eq. (3)
    
    % Wavenumbers 
    kp = -omega/(c + u);
    km =  omega/(c - u);
    % Reconstruct signals:
    ph = F*exp(1i*kp*x) + G*exp(1i*km*x);
    uh = 1/(rho*c)*( F*exp(1i*kp*x) - G*exp(1i*km*x));
end