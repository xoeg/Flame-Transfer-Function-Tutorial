function [] = fn_animate(Sim,stepsize)
    %  Function Animate the simulation
    
    % Read simulation data:
    ts = Sim.ts;
    Gs = Sim.Gs;
    Hs = Sim.Hs;
    Qs = Sim.Qs;
    Xis = Sim.Xis;
    Geom = Sim.Geom;
    Mean = Sim.Mean;
    Forcing = Sim.Forcing;
    Measurement = Sim.Measurement;
    
    % Discretize space
    Lu = Geom.Lu;
    Lb = Geom.Lb;
    b  = Geom.b;
    r  = Geom.r; 
    x  = linspace(-Lu,Lb,1000);
    rs = linspace(-b,b,2);
    [X,R] = meshgrid(x,rs);
    
    dt = ts(2) - ts(1);
    Tu  = Mean.Tu;
    Tau = Mean.Tau;
    J   = ceil((Tau + Tu)/(dt*100))*100;   
    
    cmap = cbrewer2('div','RdBu',1000);
        
    % Max: p, u, Q estimate:
    Gmax = max(abs(Gs));
    Hmax = max(abs(Hs));
    pmax = 2*max(Gmax,Hmax);
    
    U1 = Mean.u1 + 1/(Mean.rho1*Mean.c1)*(2*Gmax);
    U2 = Mean.u2 + 1/(Mean.rho2*Mean.c2)*(2*Hmax);
    umax = max(U1,U2);
    
    U1 = Mean.u1 - 1/(Mean.rho1*Mean.c1)*(2*Gmax);
    U2 = Mean.u2 - 1/(Mean.rho2*Mean.c2)*(2*Hmax);
    umin = min(U1,U2);
    qmax = 1.1*max(Qs);
    
    
    fig = figure;
    fig.Visible  = 'on';
    fig.Position = [100 100 1600*3 900*3];
    
    m = 1;
    Ts = [];
    P_mic = [];
    Us = [];
    
    frt = 24;           % Frame rate
    v1 = VideoWriter('MMM','Motion JPEG AVI');
    v1.FrameRate = frt;
    v1.Quality = 100;
    open(v1)
    
    
    % Time march animation:
    for j = 1:stepsize:length(ts)
        t = ts(j);
        % Compute pressure modeshape:
        Acoustic_only = true;
        p = zeros(size(x));
        parfor k = 1:length(x)
            [~,pk] = fn_modehshapes(t,x(k),ts,Gs,Hs,Mean,Geom,Acoustic_only);
            p(k) = pk;
        end 
        P = repmat(p,[2,1]);
        
        % Plot
        clf
        % Remove axes:
        ax = gca;
        ax.XTick = [];
        ax.YTick = [];
        ax.XLabel = [];
        ax.YLabel = [];
        ax.YColor = [1 1 1];
        ax.XColor = [1 1 1];
        box off

        ax1 = axes();
        ax1.Position = [0.027 0.13 0.98 0.26];
        
        hold on
        contourf(X,R,P,100,'Edgecolor','none')
        colormap(cmap);
        c = colorbar();
        caxis([-pmax,pmax])
        c.Location = 'southoutside';
        c.Label.String = '$p''(t)$ (Pa)';
        c.Label.Interpreter = 'Latex';
        c.TickLabelInterpreter = 'Latex';
        fn_plot_config(Forcing,Measurement)
        % Plot the flame front
        if Mean.Qbar > 0
            Xi_Tau = fn_flame_front_hist(t-Tau,Xis,ts,j,J,Mean);
            plot(Xi_Tau, r,'r','Linewidth',2)
            plot(Xi_Tau,-r,'r','Linewidth',2)
        end
        xlim([-Lu*1.07,Lb*1.01])
        ax1.YTick = [];
        ax1.YLabel = [];
        ax1.YColor = [1 1 1];
        box off
        
        ax1.FontSize = 16;
        fontsz = ax1.FontSize;
        
        % Plot measurements:
        if strcmp(Measurement.PMT,'ON')
            ax2 = axes();
            ax2.Position = [0.55 0.5 0.3 0.3];
            
            yellow = color_palette('Yl');
            plot(ts(1:j),Qs(1:j),'Color',yellow,'LineWidth',2)
            xlabel('$t$ (s)','Interpreter','Latex')
            ylabel('$Q(t)$ (W)','Interpreter','Latex')
            ax2.TickLabelInterpreter = 'Latex';
            jmin = max(1,j-1000);
            if j > 2
                xlim([ts(jmin) ts(j)])
            end
            ylim([0 qmax])
            box off
            ax2.FontSize = fontsz;
            box on
        end
        
        if strcmp(Measurement.Mic,'ON')
            
            ax3 = axes();
            ax3.Position = [0.18 0.45 0.3 0.15];
            
            Acoustic_only = true;
            x_mic = Measurement.Mic_Pos;
            N_mic = length(x_mic);
            for k = 1:N_mic
                [~,pp] = fn_modehshapes(t,x_mic(k),ts,Gs,Hs,Mean,Geom,Acoustic_only);
                P_mic(k,m) = pp;
                Ts(m) = t;
            end
            
            blus = color_palette('Bl');
            hold on
            for k = 1:N_mic
                plot(Ts,P_mic(k,:),'Color',blus*k/N_mic,'LineWidth',2)
                legendCell{k} = ['$x = ' num2str(round(x_mic(k),2)),'$ m'];
            end
            plot([0 Ts(end)*2],[0 0],'k')
            lg = legend(legendCell);
            lg.Interpreter = 'Latex';
            lg.Location = 'northoutside';
            lg.Orientation = 'horizontal';
            xlabel('$t$ (s)','Interpreter','Latex')
            ylabel('$p''(t)$ (Pa)','Interpreter','Latex')
            ax3.TickLabelInterpreter = 'Latex';
            jmin = max(1,j-1000);
            if j > 2
                xlim([ts(jmin) ts(j)])
            end
            ylim([-1,1]*pmax)
            box off
            ax3.FontSize = fontsz;
            box on
        end
        
        if strcmp(Measurement.HWA,'ON')
            x_hw = Measurement.HWA_Pos;
            Acoustic_only = false;
            up = fn_modehshapes(t,x_hw,ts,Gs,Hs,Mean,Geom,Acoustic_only);
            Us(m) = up;
            
            ax4 = axes();
            ax4.Position = [0.18 0.70 0.3 0.15];
            
            reds = color_palette('Rd');
            hold on
            plot(Ts,Us,'Color',reds,'LineWidth',2)
            plot([0 Ts(end)*2],[0 0],'k')
            xlabel('$t$ (s)','Interpreter','Latex')
            ylabel('$u(t)$ (m/s)','Interpreter','Latex')
            ax4.TickLabelInterpreter = 'Latex';
            jmin = max(1,j-1000);
            if j > 2
                xlim([ts(jmin) ts(j)])
            end
            ylim([umin, umax])
            box off
            ax4.FontSize = fontsz;
            box on
            
        end
        m = m + 1;
        drawnow
        pause(0.0001)
        
        % Frame
        Mv = getframe(gcf);
        writeVideo(v1,Mv)
        
        j/length(ts)
        
    
    end
    close(v1);

end


