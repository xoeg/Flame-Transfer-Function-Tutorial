%% Modeshapes
function [up,pp,rhop] = fn_modehshapes(t,x,ts,Gs,Hs,Mean,Geom,Acoustic_only)
    % Computes the modeshapes at t, x
    
    if x <= 0
        % Retrieve variables:
        A1 = Geom.A1;
        M1 = Mean.M1;
        c1 = Mean.c1;
        u1 = Mean.u1;
        p1 = Mean.p1;
        rho1 = Mean.rho1;
        Tu = Mean.Tu;
        % Definition of f and wave functions
        if Geom.Speaker
            Vu     = Mean.Force.Vu;
            omega_u = Mean.Force.omega_u;
            fv = @(t) Vu*c1/(A1*(1 + M1))*cos(omega_u*(t - Tu/2));
        else
            fv = @(t) zeros(size(t));
        end
        g = @(t) fv_interp(t,ts,Gs);
        f = @(t) (1 - M1)/(1 + M1) * g(t - Tu) + fv(t);
        F = f(t - x/(c1 + u1));
        G = g(t + x/(c1 - u1));
        % Modeshapes
        if Acoustic_only
            pp   = F + G;
            up   = 1/(rho1*c1)*(F - G);
            rhop = rho1 + 1/(c1^2) * (F + G);
        else
            pp   = p1 + F + G;
            up   = u1 + 1/(rho1*c1)*(F - G);
            rhop = rho1 + 1/(c1^2) * (F + G);
        end
    else
        % Retrieve variables:
        c2 = Mean.c2;
        u2 = Mean.u2;
        p2 = Mean.p2;
        rho2 = Mean.rho2;
        Td = Mean.Td;
        % Definition of j and wave functions
        if Geom.Speaker_Downstream
            Vd    = Mean.Force.Vd;
            omega_d = Mean.Force.omega_d;
            fd = @(t) Vd*cos(omega_d*(t - Td/2));
        else
            fd = @(t) zeros(size(t));
        end
        h = @(t) fv_interp(t,ts,Hs);
        j = @(t) - h(t - Td) + fd(t);
        H = h(t - x/(c2 + u2));
        J = j(t + x/(c2 - u2));
        % Modeshapes
        if Acoustic_only
            pp = H + J;
            up = 1/(rho2*c2)*(H - J);
            rhop = NaN; % Missing entropy waves
        else
            pp = p2 + H + J;
            up = u2 + 1/(rho2*c2)*(H - J);
            rhop = NaN; % Missing entropy waves
        end
    end
end

%% Function definitions:
function v = fv_interp(t,ts,Vs)
    % Interpolate
    if t < 0
        v = 0;
    else
        v =  interp1(ts,Vs,t);
    end
end