%% Generate Animations:
close all; clear all; clc

% Define inlet and outlet conditions. 
u1 = 20;       % Inlet velocity (m/s)
T1 = 293;      % Inlet temperatre (m/s)
p2 = 101300;   % Outlet pressure (Pa)
Qbar = 0; % Heat release rate (W)

% Define Forcing Conditions:
Forcing.Speaker = 'ON';
Forcing.Voltage = 1; % (V)
Forcing.Frequency = 200; % (Hz)
% Measurement probes:

% Sampling frequency: 
Measurement.Fs = 5000; % (Hz)
% Photomultiplier:
Measurement.PMT = 'OFF';
% Differential Microphones:
Measurement.Mic = 'ON';
Measurement.Mic_Pos  = [-0.2,-0.4]; % (m) relative to the configuration
% Hot wire anemometry:
Measurement.HWA = 'ON';
Measurement.HWA_Pos = 0.2; % (m) 

% Length of the measurements:
Measurement.Time = 2; % (s)

% Preview set up:
% figure('Visible','on')
% fn_plot_config(Forcing,Measurement)

% Simulation
[Mean,Geom,Signals,Sim] = Run_Simulation(u1,T1,p2,Qbar,Forcing,Measurement);

% hold on
% plot(Sim.ts,Sim.Gs)
% plot(Sim.ts,Sim.Hs)
% return

%% Signal reconstruction:
% Multi-Microphone-Method:
[F,G,omega] = fn_MMM_comp(Mean,Signals,Measurement);
% 1) Compute the values of uh1 and ph1:
[ph1,uh1] = fn_acoustic_field_comp(Mean,F,G,omega,0);
% 2) Compute F and G after the area change:
[F2,G2] = fn_Area_Change_comp(Mean,Geom,ph1,uh1,omega,0);
% 3) Propagate the signal:
[~,uh] = fn_acoustic_field_comp(Mean,F2,G2,omega,0.2);
% Retrieve the mean velocity
t = Sim.ts;
% Time domain (harmonic) signal
ur = Mean.u2 + real(uh*exp(1i*omega*t));
% Save and pass
Sim.Measurement.REC = 'ON';
Sim.ur = ur;


%% Animate:
stepsize = 1;
% fn_animate(Sim,stepsize)
MAX.umax = 17;
MAX.umin = 13;
MAX.pmax = 2000;
tstart = 1;
tend = 1.005;
fn_gifgen(Sim,stepsize,tstart,tend,MAX)
 
