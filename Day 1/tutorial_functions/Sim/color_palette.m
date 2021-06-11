function [clr] = color_palette(color)
    % Color palete

    Bl = [0,102,204]/255;
    Gr = [0.235294117647059,0.682352941176471,0.63921568627451];
    Yl = [0.964705882352941,0.835294117647059,0.33725490196078];
    Rd = [0.929411764705882,0.333333333333333,0.23137254901960];
    
    switch color
        case 'Bl'; clr = Bl;
        case 'Gr'; clr = Gr;
        case 'Yl'; clr = Yl;
        case 'Rd'; clr = Rd;
    end
end
