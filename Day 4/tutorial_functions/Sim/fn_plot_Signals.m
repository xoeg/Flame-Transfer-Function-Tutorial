function [] = fn_plot_Signals(Signals,Measurement,Zoom_in_time_steps)
    % Function to plot all the signals
    
    % Determine the number of plots
    N = 0;
    fn = fieldnames(Signals);
    for k = 1:numel(fn)
        if strcmp(fn{k},'Forcing')
            N = N + size(Signals.Forcing,1);
        elseif strcmp(fn{k},'Forcing_Downstream')
            N = N + size(Signals.Forcing_Downstream,1);
        elseif strcmp(fn{k},'P_mic')
            N = N + size(Signals.P_mic,1);
        elseif strcmp(fn{k},'U_hw')
            N = N + size(Signals.U_hw,1);
        elseif strcmp(fn{k},'PMT')
            N = N + size(Signals.PMT,1);
        end
    end
    % Zoom in region:
    t = Signals.t;
    nn = length(t);
    XL = [t(round(nn/2)) t(round(nn/2) + Zoom_in_time_steps)];
    % Plot the figure
    figure('Color','w')
    jj = 1;
    for k = 1:numel(fn)
        if strcmp(fn{k},'Forcing')
            Nm = size(Signals.Forcing,1);
            X = Signals.Forcing;
            for j = 1:Nm
                clr = 'k';
                subplot(N,3,[jj,jj+1])
                plot(t,X(j,:),'Color',clr)
                xlabel('$t$ (s)','Interpreter','Latex')
                ylabel('$F(t)$ (V)','Interpreter','Latex')
                ax = gca;
                ax.TickLabelInterpreter = 'Latex';
                box on
                % Zoom in
                subplot(N,3,jj+2)
                plot(t,X(j,:),'Color',clr)
                xlabel('$t$ (s)','Interpreter','Latex')
                ylabel('$F(t)$ (V)','Interpreter','Latex')
                ax = gca;
                ax.TickLabelInterpreter = 'Latex';
                box on
                xlim(XL)
                jj = jj + 3;
            end
        elseif strcmp(fn{k},'Forcing_Downstream')
            Nm = size(Signals.Forcing_Downstream,1);
            X = Signals.Forcing_Downstream;
            for j = 1:Nm
                clr = 'k';
                subplot(N,3,[jj,jj+1])
                plot(t,X(j,:),'Color',clr)
                xlabel('$t$ (s)','Interpreter','Latex')
                ylabel('$F(t)$ (V)','Interpreter','Latex')
                ax = gca;
                ax.TickLabelInterpreter = 'Latex';
                box on
                % Zoom in
                subplot(N,3,jj+2)
                plot(t,X(j,:),'Color',clr)
                xlabel('$t$ (s)','Interpreter','Latex')
                ylabel('$F(t)$ (V)','Interpreter','Latex')
                ax = gca;
                ax.TickLabelInterpreter = 'Latex';
                box on
                xlim(XL)
                jj = jj + 3;
            end
        elseif strcmp(fn{k},'P_mic')
            Nm = size(Signals.P_mic,1);
            t = Signals.t;
            X = Signals.P_mic;
            for j = 1:Nm
                clr = color_palette('Bl');
                subplot(N,3,[jj,jj+1])
                plot(t,X(j,:),'Color',clr)
                xlabel('$t$ (s)','Interpreter','Latex')
                ylabel('$p''(t)$ (Pa)','Interpreter','Latex')
                ax = gca;
                ax.TickLabelInterpreter = 'Latex';
                box on
                x = Measurement.Mic_Pos(j);
                legend(['$x = ',num2str(x),'$ m'],'Interpreter','Latex')
                % Zoom in
                subplot(N,3,jj+2)
                plot(t,X(j,:),'Color',clr)
                xlabel('$t$ (s)','Interpreter','Latex')
                ylabel('$p''(t)$ (Pa)','Interpreter','Latex')
                ax = gca;
                ax.TickLabelInterpreter = 'Latex';
                box on
                xlim(XL)
                jj = jj + 3;
            end
        elseif strcmp(fn{k},'U_hw')
            clr = color_palette('Rd');
            X = Signals.U_hw;
            subplot(N,3,[jj,jj+1])
            plot(t,X,'Color',clr)
            xlabel('$t$ (s)','Interpreter','Latex')
            ylabel('$u(t)$ (m/s)','Interpreter','Latex')
            ax = gca;
            ax.TickLabelInterpreter = 'Latex';
            box on
            x = Measurement.HWA_Pos;
            legend(['$x = ',num2str(x),'$ m'],'Interpreter','Latex')
            % Zoom in 
            subplot(N,3,jj+2)
            plot(t,X,'Color',clr)
            xlabel('$t$ (s)','Interpreter','Latex')
            ylabel('$u(t)$ (m/s)','Interpreter','Latex')
            ax = gca;
            ax.TickLabelInterpreter = 'Latex';
            box on
            xlim(XL)
            jj = jj + 3;
        elseif strcmp(fn{k},'PMT')
            clr = color_palette('Yl');
            X = Signals.PMT;
            subplot(N,3,[jj,jj+1])
            plot(t,X,'Color',clr)
            xlabel('$t$ (s)','Interpreter','Latex')
            ylabel('$I(t)$ ','Interpreter','Latex')
            ax = gca;
            ax.TickLabelInterpreter = 'Latex';
            box on
            % Zoom in 
            subplot(N,3,jj+2)
            plot(t,X,'Color',clr)
            xlabel('$t$ (s)','Interpreter','Latex')
            ylabel('$I(t)$ ','Interpreter','Latex')
            ax = gca;
            ax.TickLabelInterpreter = 'Latex';
            box on
            xlim(XL)
            jj = jj + 3;
        end
    end
end