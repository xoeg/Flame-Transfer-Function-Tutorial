function [] = fn_plot_config(Forcing,Measurement)
    % Function that plots the geometric configuration:

    % Set up =============================================================
    % Set up measurements:
    % Hotwire:
    if strcmp(Measurement.HWA,'ON') && length(Measurement.HWA_Pos) > 1
        error('Our Lab has only one hot wire anemometer.')
    end
    % Load the geometry
    order = 2;
    Geom = fn_define_geometry(order);
    % Set up Forcing:
    if strcmp(Forcing.Speaker,'ON')
        Geom.Speaker = true;
    else
        Geom.Speaker = false;
    end
    if strcmp(Forcing.Speaker_Downstream,'ON')
        Geom.Speaker_Downstream = true;
    else
        Geom.Speaker_Downstream = false;
    end
    % Plot ===============================================================
    Lu = Geom.Lu;
    Lb = Geom.Lb;
    b  = Geom.b;
    a  = Geom.a;
    h  = b/5;       % Thickness of duct:
    % Plot the ducts
    hold on   
    patch([-Lu,Lb,Lb,-Lu,-Lu],[b,b,b+h,b+h,b],[0.8 0.8 0.8])
    patch([-Lu,Lb,Lb,-Lu,-Lu],[-b,-b,-b-h,-b-h,-b],[0.8 0.8 0.8])
    patch([-Lu,0,0,-Lu,-Lu],[-a,-a,a,a,-a],[0.8 0.8 0.8])
    axis equal
    s = (Lu + Lb) *0.05;
    xlim([-Lu-s,Lb+s])
    ylim([-b-s,b+s])
    % Plot Speaker:
    if Geom.Speaker
        h = b;
        r = [b,-b,-b/2,-b/2,b/2,b/2,b];
        x = [2*h,2*h,h,0,0,h,2*h];
        
        patch(x-Lu-2*h,r,[0.2 0.2 0.2])
        
        xlim([-Lu-3*h-s,Lb+s])
    end
    % Plot Downstream speaker:
    if Geom.Speaker_Downstream
        h = b;
        r = [b,-b,-b/2,-b/2,b/2,b/2,b];
        x = [2*h,2*h,h,0,0,h,2*h];
        
        patch(Lb-x+2*h,r,[0.2 0.2 0.2])
        
        xlim([-Lu-3*h-s,Lb+s+2*h])
    end
    % Plot the PMT
    if strcmp(Measurement.PMT,'ON')
        clr = color_palette('Yl');
        h = b/2;
        x0 = 0.075-h/2;
        y0 = b + 2*b/5+h;
        patch(x0+[0 2*h 2*h 0 0],y0+[0 0 3*h 3*h 0],clr)
        x0 = 0.075;
        y0 = b + 2*b/5;
        patch(x0+[0 h h 0 0],y0+[0 0 h h 0],[0.3 0.3 0.3])
    end
    % Plot the Microphones:
    if strcmp(Measurement.Mic,'ON')
        x_mic = Measurement.Mic_Pos;
        clr = color_palette('Bl');
        h = b/5;
        for j = 1:length(x_mic)
            x0 = x_mic(j)-h/2;
            y0 = b;
            patch(x0+[0 h h 0 0],y0+[0 0 3*h 3*h 0],clr)
        end
    end
     % Plot the HW:
    if strcmp(Measurement.HWA,'ON')
        x_hw = Measurement.HWA_Pos;
        clr = color_palette('Rd');
        h = b/5;
        x0 = x_hw-h/2;
        y0 = b;
        patch(x0+[0 h h 0 0],y0+[0 0 3*h 3*h 0],clr)
    end
    
    
    % Format
    box on
    xlabel('$x$ (m)','Interpreter','Latex')
    ylabel('$r$ (m)','Interpreter','Latex')
    ax = gca;
    ax.TickLabelInterpreter = 'Latex';
    
    
end