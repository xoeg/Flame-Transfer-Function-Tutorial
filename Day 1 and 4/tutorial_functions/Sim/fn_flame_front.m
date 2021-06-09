function [dxi_dt]  = fn_flame_front(t,xi,ts,Xis,Gs,Hs,j,J,Mean,Geom)
    % Computes the Flame Front as in Eq. (2.8)
    
    % Retreive
    Su = Mean.Su;
    % Wave functions
    [gt,~,gt_Tu] = fun_wavefun(t,ts,Xis,Gs,Hs,j,J,Mean,Geom);
    % Velocity definitions
    Ut = fn_flame_holder_vel(t,gt,gt_Tu,Mean,Geom);
    dxi_dr = fn_dxi_dr(xi,Ut,Mean,Geom);
    % Derivative of flame front wrt t:
    dxi_dt = Ut*ones(size(xi)) - Su*sqrt(1 + (dxi_dr).^2);
end