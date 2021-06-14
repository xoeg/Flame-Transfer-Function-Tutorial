%% Wrapper function:
close all
clear all
clc

WS = load('FTF_8kW.mat');
FTF = WS.FTF_8kw;
f   = WS.f_8kw;

% Oder the FTF
[f,ix] = unique(f);
FTF = FTF(ix);
omegas = 2*pi*f;


FB = [2000 9000]; %Frequency Band of modulations

flag_plot = true;
[V] = fn_Fit_Distr(omegas,FTF,FB,flag_plot);

%%
close all
clear all
clc

WS = load('FTF_12kW.mat');
FTF = WS.FTF_12kw;
f   = WS.f_12kw;

% Oder the FTF
[f,ix] = unique(f);
FTF = FTF(ix);
omegas = 2*pi*f;


FB = [2000 12000]; %Frequency Band of modulations

flag_plot = true;
[V] = fn_Fit_Distr(omegas,FTF,FB,flag_plot);
