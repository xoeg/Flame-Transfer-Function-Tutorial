%----------- Amplitude dynamics and statistics of Van der Pol oscillator -----------------
%----------- Nicolas Noiray. last update: nov. 2016
% Matlab script and simulink model used to generate the figure 4
% of the paper entitled "Linear Growth Rate Estimation From Dynamics 
% and Statistics of Acoustic Signal Envelope in Turbulent Combustors" 
% from Nicolas Noiray and published in 2017 in the 
% Journal of Engineering for Gas Turbines and Power


clear all
close all

% Parameters of the stochastic Van der Pol oscillator
omega=2*pi*120;
nu=4;
kappa=3.4;
Gamma=1e7;
% Theoretical PDF
a_vec=[0.01:0.01:0.1 0.2:0.1:12]; % Amplitude vector
pdf_A=a_vec.* exp(4*omega^2/Gamma *(nu*a_vec.^2/2-kappa*a_vec.^4/32));
pdf_A_norm=pdf_A./sum(pdf_A.*[0 diff(a_vec)]); % Normalized PDF

figure
plot(a_vec,pdf_A)
xlim([0 6])
xlabel('A','fontsize',16)
ylabel('PDF(A)','fontsize',16)


%% 
% Compute stochastic bifurcation diagram 
nu=-15:0.1:15;
for ii=1:length(nu);
    pdf_Ab(ii,:)=a_vec.* exp(4*omega^2/Gamma *(nu(ii)/2*a_vec.^2-kappa/32*a_vec.^4)) ...
    /sum(a_vec.* exp(4*omega^2/Gamma *(nu(ii)/2*a_vec.^2-kappa/32*a_vec.^4)).*[0 diff(a_vec)]);
    pot_Ab(ii,:)=-(nu(ii)/2*a_vec.^2-kappa/32*a_vec.^4+Gamma/4/omega^2*log(a_vec));
end
[Amp,Nu]=meshgrid(a_vec,nu);

%%
% Perform numerical simulations of the Van der Pol oscillator forced by
% white noise using Simulink for 3 values of the linear growth rate:
% linearly stable (negative growth rate)
% marginaly stable (growth rate=0)
% linearly unstable (positive growth rate)

alpha_i=[48 60 72];
for ii=1:3
T=1*60; % Duration of the simulation (s)
Fs=1000; % Sampling frequency (Hz)
wn=2*pi*120; % oscillator frequency
alpha=alpha_i(ii);% linear damping
beta=60; % linear gain
nu0=(beta-alpha)/2; % linear growth/decay rate  (positive/negative)
random_seed=floor(1e3*rand(1)); % Seed source block
sim VDP_OSCILLATOR
ind_trans=10000;  % index to remove transient part of signal before computing the PDF from the histogram of the simulated data
xi=simout.signals.values(ind_trans:end,1);  % 
eta=simout.signals.values(ind_trans:end,2);
time=simout.time(1:end-ind_trans+1);
nsamp=length(time);
Amp0=abs(hilbert(eta));  % extraction of the enveloppe of the signal using the Hilbert transform (other enveloppe extraction methods exist)
[nya,cya]=hist(Amp0,100); % compute histogram of enveloppe
nya=nya/trapz(cya,nya); % deduce PDF by normalizing the histogram
pdf_tt{ii}.pdf=nya;
pdf_tt{ii}.a=cya;
end

%%
% Plot the results

figure('position',[16   424   865   276]);
subplot('position',[0.05    0.2    0.2    0.7]);
cm=hot(256);cm=cm(end:-1:1,:);
contourf(Nu,Amp,pdf_Ab,'linecolor','none');colormap(cm);
hold on
st_point=sqrt(8*nu/kappa).*(nu>0);
plot(nu,st_point,'k','linewidth',4)
plot(-6*[1 1],[0 12],'--k')
plot(0*[1 1],[0 12],'--k')
plot(6*[1 1],[0 12],'--k')
xlabel('\nu (rad/s)','fontsize',16)
ylabel('A','fontsize',16)
caxis([0 1.2]);%colorbar
set(gca,'fontsize',16,...
    'box','on')
ylim([0 12])

subplot('position',[0.36    0.46    0.2    0.45])
hold on
set(gca,'ytick',[-20 0 20],'xtick',[0 2 4 6],'xticklabel','','fontsize',16,...
    'box','on')
plot(a_vec,-nu(find(nu==-6))/2*a_vec.^2,'--','color',[0.7 0 0],'linewidth',2)
plot(a_vec,+kappa/32*a_vec.^4,'--','color',[0 0 0.7],'linewidth',2)
plot(a_vec,-Gamma/4/omega^2*log(a_vec),'--','color',[0 0.7 0],'linewidth',2)
plot(a_vec,pot_Ab(nu==-6,:),'k','linewidth',3)
ylabel('Potential','fontsize',16)
title('\nu=-6 rad/s','fontsize',16)
xlim([0 7]);ylim([-30 30]);grid on

subplot('position',[0.57    0.46    0.2    0.45])
hold on
set(gca,'ytick',[-20 0 20],'xtick',[0 2 4 6],'yticklabel','','xticklabel','','fontsize',16,...
    'box','on')
plot(a_vec,-nu(find(nu==0))/2*a_vec.^2,'--','color',[0.7 0 0],'linewidth',2)
plot(a_vec,+kappa/32*a_vec.^4,'--','color',[0 0 0.7],'linewidth',2)
plot(a_vec,-Gamma/4/omega^2*log(a_vec),'--','color',[0 0.7 0],'linewidth',2)
plot(a_vec,pot_Ab(nu==0,:),'k','linewidth',3)
ylabel('')
title('\nu=0 rad/s','fontsize',16)
xlim([0 7]);ylim([-30 30]);grid on

subplot('position',[0.78    0.46    0.2    0.45])
hold on
set(gca,'ytick',[-20 0 20],'xtick',[0 2 4 6],'yticklabel','','xticklabel','','fontsize',16,...
    'box','on')
plot(a_vec,-nu(find(nu==6))/2*a_vec.^2,'--','color',[0.7 0 0],'linewidth',2)
plot(a_vec,+kappa/32*a_vec.^4,'--','color',[0 0 0.7],'linewidth',2)
plot(a_vec,-Gamma/4/omega^2*log(a_vec),'--','color',[0 0.7 0],'linewidth',2)
plot(a_vec,pot_Ab(nu==6,:),'k','linewidth',3)
xlabel('');ylabel('')
title('\nu=6 rad/s','fontsize',16)
xlim([0 7]);ylim([-30 30]);grid on


subplot('position',[0.36    0.2   0.2    0.22])
hold on
patch([a_vec a_vec],...
    [zeros(size(a_vec)) pdf_Ab(nu==-6,:)],0.8*ones(1,3))  % theoretical PDF
plot(pdf_tt{3}.a,pdf_tt{3}.pdf,'b','linewidth',2)       % PDF from simulated data
xlabel('A','fontsize',16)
ylabel('PDF','fontsize',16)
set(gca,'fontsize',16,'xtick',[0 2 4 6],...
    'box','on')
xlim([0 7]);ylim([0 1]);grid on

subplot('position',[0.57    0.2   0.2    0.22])
hold on
patch([a_vec a_vec],...
    [zeros(size(a_vec)) pdf_Ab(nu==0,:)],0.8*ones(1,3))  % theoretical PDF
plot(pdf_tt{2}.a,pdf_tt{2}.pdf,'b','linewidth',2)  % PDF from simulated data
xlabel('A','fontsize',16)
set(gca,'yticklabel','','fontsize',16,'xtick',[0 2 4 6],...
    'box','on')
xlim([0 7]);ylim([0 1]);grid on

subplot('position',[0.78    0.2   0.2    0.22])
hold on
patch([a_vec a_vec],...
    [zeros(size(a_vec)) pdf_Ab(nu==6,:)],0.8*ones(1,3))  % theoretical PDF
plot(pdf_tt{1}.a,pdf_tt{1}.pdf,'b','linewidth',2)  % PDF from simulated data
xlabel('A','fontsize',16)
set(gca,'yticklabel','','fontsize',16,'xtick',[0 2 4 6],...
    'box','on')
xlim([0 7]);ylim([0 1]);grid on




