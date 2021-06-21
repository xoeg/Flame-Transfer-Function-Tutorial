function [freqs,X22] = fn_plot_BTM(fs, X11e, X12e, X21e, X22e, str)
    % Function that plots the BTM 
    if strcmp(str,'BTM')
        WS = load('BTM_1.mat');
        freqs = WS.freqs;
        X11 = WS.B11;
        X12 = WS.B12;
        X21 = WS.B21;
        X22 = WS.B22;
    elseif strcmp(str,'FTM_1') || strcmp(str,'FTM_2')
        WS = load('BTM_1.mat');
        B11 = WS.B11;
        B12 = WS.B12;
        B21 = WS.B21;
        B22 = WS.B22;
        WS = load([str,'.mat']);
        freqs = WS.freqs;
        T11 = WS.T11;
        T12 = WS.T12;
        T21 = WS.T21;
        T22 = WS.T22;
        [X11,X12,X21,X22] = deal(zeros(size(T11)));
        for j = 1:length(T11)
            B = [B11(j) B12(j);B21(j) B22(j)];
            T = [T11(j) T12(j);T21(j) T22(j)];
            F = B\T;
            X11(j) = F(1,1);
            X12(j) = F(1,2);
            X21(j) = F(2,1);
            X22(j) = F(2,2);
        end
    end
    % Remove the values of X at low frequencies to avoid the spike:
    freqs = freqs(2:end);
    X11 = X11(2:end);
    X12 = X12(2:end);
    X21 = X21(2:end);
    X22 = X22(2:end);
    
    M = X11;
    Me = X11e;
    subplot(4,2,1)
    plot(freqs,abs(M))
    hold on
    plot(fs,abs(Me),'o')
    if strcmp(str,'BTM')
        ylabel('$|B_{11}|$','Interpreter','Latex')
    else
        ylabel('$|F_{11}|$','Interpreter','Latex')
    end
    box on
    grid on
    ax = gca;
    ax.XTickLabels = [];
    ylim([0 1.1*max(abs(M))])
    ax.TickLabelInterpreter = 'Latex';
    
    subplot(4,2,3)
    plot(freqs,angle(M))
    hold on
    plot(fs,angle(Me),'o')
    xlabel('$f$ [Hz]','Interpreter','Latex')
    if strcmp(str,'BTM')
        ylabel('$\angle(B_{11})$','Interpreter','Latex')
    else
        ylabel('$\angle(F_{11})$','Interpreter','Latex')
    end
    box on
    grid on
    ax = gca;
    ax.TickLabelInterpreter = 'Latex';
    

    M = X12;
    Me = X12e;
    subplot(4,2,2)
    plot(freqs,abs(M))
    hold on
    plot(fs,abs(Me),'o')
    if strcmp(str,'BTM')
        ylabel('$|B_{12}|$','Interpreter','Latex')
    else
        ylabel('$|F_{12}|$','Interpreter','Latex')
    end
    box on
    grid on
    ax = gca;
    ax.XTickLabels = [];
    ylim([0 1.1*max(abs(M))])
    ax.TickLabelInterpreter = 'Latex';
    
    subplot(4,2,4)
    plot(freqs,angle(M))
    hold on
    plot(fs,angle(Me),'o')
    xlabel('$f$ [Hz]','Interpreter','Latex')
    if strcmp(str,'BTM')
        ylabel('$\angle(B_{12})$','Interpreter','Latex')
    else
        ylabel('$\angle(F_{12})$','Interpreter','Latex')
    end
    box on
    grid on
    ax = gca;
    ax.TickLabelInterpreter = 'Latex';
    
    M = X21;
    Me = X21e;
    subplot(4,2,5)
    plot(freqs,abs(M))
    hold on
    plot(fs,abs(Me),'o')
    if strcmp(str,'BTM')
        ylabel('$|B_{21}|$','Interpreter','Latex')
    else
        ylabel('$|F_{21}|$','Interpreter','Latex')
    end
    box on
    grid on
    ax = gca;
    ax.XTickLabels = [];
    ylim([0 1.1*max(abs(M))])
    ax.TickLabelInterpreter = 'Latex';
    
    subplot(4,2,7)
    plot(freqs,angle(M))
    hold on
    plot(fs,angle(Me),'o')
    xlabel('$f$ [Hz]','Interpreter','Latex')
    if strcmp(str,'BTM')
        ylabel('$\angle(B_{21})$','Interpreter','Latex')
    else
        ylabel('$\angle(F_{21})$','Interpreter','Latex')
    end
    box on
    grid on
    ax = gca;
    ax.TickLabelInterpreter = 'Latex';

    M = X22;
    Me = X22e;
    subplot(4,2,6)
    plot(freqs,abs(M))
    hold on
    plot(fs,abs(Me),'o')
    if strcmp(str,'BTM')
        ylabel('$|B_{22}|$','Interpreter','Latex')
    else
        ylabel('$|F_{22}|$','Interpreter','Latex')
    end
    box on
    grid on
    ax = gca;
    ax.XTickLabels = [];
    ylim([0 1.1*max(abs(M))])
    ax.TickLabelInterpreter = 'Latex';
    
    subplot(4,2,8)
    plot(freqs,angle(M))
    hold on
    plot(fs,angle(Me),'o')
    xlabel('$f$ [Hz]','Interpreter','Latex')
    if strcmp(str,'BTM')
        ylabel('$\angle(B_{22})$','Interpreter','Latex')
    else
        ylabel('$\angle(F_{22})$','Interpreter','Latex')
    end
    box on
    grid on
    ax = gca;
    ax.TickLabelInterpreter = 'Latex';

end