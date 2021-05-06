function [F,G] = fn_Area_Change(Mean,Geom,ph1,uh1,omega,x0)
    % Function that copmutes the values of F and G after an area
    % change. Receives as an input the "Mean" and "Geom" structures, the
    % values of ph1 and uh1, the angular frequency "omega" and the position
    % "x0". It outputs the amplitude values of "F" and "G" after the area
    % change.
    
    % Retrieve variables:
    A1 = Geom.A1;
    A2 = Geom.A2;
    rho = Mean.rho2;
    c   = Mean.c2;
    u   = Mean.u2;
    % ====================================================================
    % Acoustic field after the jump --------------------------------------
    % Compute the values of ph2 and uh2 Using Eq. (5)

    
    % Compute the values of F and G using Eq. (3) ------------------------

    
    
end