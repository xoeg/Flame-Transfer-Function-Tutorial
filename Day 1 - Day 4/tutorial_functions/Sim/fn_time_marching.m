function [ts,Gs,Hs,Xis,Qs] = fn_time_marching(t_end,dt,Geom,Mean)
    % Function to time march a simulation:
    
    % Retreive variables:
    Tau = Mean.Tau;
    Tu = Mean.Tu;
    N = Geom.N;
    
    % Time marching paramters:
    Nt = ceil(t_end/dt);
    Nt = Nt + 1;
    J  = ceil((Tau + Tu)/dt)+10;    
    
    % Initial conditions:
    xi = Mean.xibar;
    t = 0;

    % Prelocate
    Gs = zeros(Nt,1);
    Hs = zeros(Nt,1);
    Qs = zeros(Nt,1);
    ts = zeros(Nt,1);
    Xis = zeros(N,Nt);
    for j = 1:Nt
        if j == 1
            fprintf('Progress: 0%%,')
        elseif mod(j,1000) == 0
            if round(j/Nt,2)*100 < 100
                fprintf([' ',num2str(round(j/Nt,2)*100),'%%, '])
            end
        end
        % Compute g(t) and h(t)
        [gt,ht,~,Qt] = fun_wavefun(t,ts,Xis,Gs,Hs,j,J,Mean,Geom);
        Gs(j) = gt;
        Hs(j) = ht;
        Qs(j) = Qt;
        ts(j) = t;
        Xis(:,j) = xi;
        if Mean.Qbar == 0
            t = t + dt;
        else
            % Compute new xi:
            fun = @(t,xi) fn_flame_front(t,xi,ts,Xis,Gs,Hs,j,J,Mean,Geom);
            % Time march
            [xi,t] = fn_4thRunge_Kutta(t,xi,dt,fun);
            if max(abs(xi)) > 1e4
                if Geom.Speaker
                    error(['Solution not found. ',...
                           'Change operating conditions or forcing.'])
                else
                    error(['Solution not found.',...
                           'Change operating conditions.'])
                end
            end
        end
    end
end