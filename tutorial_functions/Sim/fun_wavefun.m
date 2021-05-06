function [gt,ht,gt_Tu,Qt] = fun_wavefun(t,ts,Xis,Gs,Hs,j,J,Mean,Geom)
    % Function that computes g(t) and h(t) via Eq (3.9)
    
    % Retrieve variables:
    c1 = Mean.c1;
    % Network paramters:
    Tu = Mean.Tu;
    Td = Mean.Td;
    InvX  = Mean.IX;
    InvXY = Mean.IXY;
    % Flame parameters:
    Abar = Mean.Abar;
    Qbar = Mean.Qbar;
    Tau  = Mean.Tau;
    % Geometric parameters
    A1 = Geom.A1;
    r  = Geom.r;
    % Retrieve wave amplitudes at previous states:
    gt_Tau    = fwavfun_hist(t-Tau,Gs,ts,j,J);
    gt_Tu_Tau = fwavfun_hist(t-Tu-Tau,Gs,ts,j,J);
    gt_Tu     = fwavfun_hist(t-Tu,Gs,ts,j,J);
    ht_Td     = fwavfun_hist(t-Td,Hs,ts,j,J);
    % Compute unsteady heat release:
    if Qbar == 0
        Qt = 0;
    else
        % Flame position at previous states:
        xit_Tau = fn_flame_front_hist(t-Tau,Xis,ts,j,J,Mean);
        % Velocity U(t - Tau)
        Ut_Tau = fn_flame_holder_vel(t,gt_Tau,gt_Tu_Tau,Mean,Geom);
        % dxi/dr(t - Tau)
        dxi_dr = fn_dxi_dr(xit_Tau,Ut_Tau,Mean,Geom);
        % Computes flame Area: A(t - Tau) Eq. (2.9)
        A_t_Tau = trapz(r,2*pi*r.*sqrt(1 + dxi_dr.^2));
        % Compute Heat release Q(t). Eq. (3.1)
        Qt = Qbar*A_t_Tau/Abar;
    end
    % Computes g(t) and h(t) via Eq. (3.9), where
    % X*w(t) = Y*w(t-Tau) + Fr(t), w = [g(t);h(t)]
    if Geom.Speaker
        F     = Mean.F;
        V     = Mean.Force.V;
        omega = Mean.Force.omega;
        Fv = F*V*cos(omega*(t - Tu/2));
    else
        Fv = [0;0];
    end
    Fq = [0;(Qt - Qbar)/(A1*c1)];
    Fr = Fq + Fv;
    w = InvXY * [gt_Tu;ht_Td] + InvX * Fr;
    % Output:
    gt = w(1);
    ht = w(2);
end