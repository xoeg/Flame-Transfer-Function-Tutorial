%%%%%%%%%%%%%%%%%%%% ANNULIGhT Workshop June 2021 %%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script aims at giving a better understanding of flame modeling with 
% distributed-time-delay models. Various Unit-Impulse-Response are available 
% to play with. It helps understanding the impact of the UIR coefficients on
% the FTF.
%
% Authors: Wolfgang Polifke, Malte Merk and Guillaume Fournier
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

%% Set up colors for plotting
C = colormap('lines');
close(gcf);
myBlue=C(1,:);myRed = C(2,:);myYellow=C(3,:);
myPurple=C(4,:);myGreen=C(5,:);myCyan=C(6,:);
myBlack = [0 0 0];

%% Model parameters and frequency range -> modifiable
% 
% 5 DTD are available. Change the value of the variable "DTD_Case" (from 1 to 5) to change the
% data to load and explore various models. The h-array is loaded from the .mat files. It contains 
% the coefficients of the UIR.

DTD_Case = 1;                                 % study different cases. To be changed from 1 to 5

load(['UIR_case' num2str(DTD_Case) '.mat']);  % load the UIR coefficients
name = ['Fig' num2str(DTD_Case)];             % name for output pdf

omega = (0:1:1500)';                          % frequency range of interest
dt=0.002;                                     % sampling time of the UIR
tau = dt*(0:1:length(h)-1)';                  % time vector of the UIR
tau_R = tau(end);                             % restoration time, i.e. time duration of the UIR. Only used for non-dimensionilizing 

% The variable "ind" represent the index at which the markers should appear in the phasor plot and
% the bode plot. "ind" can also be an array of multiple elements. Choose
% one or multiple indices in the omega vector.
ind = [omega(round(length(omega)/4)) omega(round(length(omega)/2))];    % for example, markers at 1/4 and 1/2 of the max frequency of interest

%% Plot definitions -> modifiable
% In the following all the "design" options for the plots are specified
%% Defintions bar plot  -> modifiable
bar_c = myBlack;                    % color of bars
min_b = -0.5; max_b = 1.5;          % y-limits of bar plot
xlabel_b = '$\tau/\tau_R$';         % xlabel in latex syntax
ylabel_b = '$h_k$';                 % ylable in latex syntax
bar_lw = 0.2;                       % width of the bars
xtick_size_b = 0.1;                 % distance between x-ticks
ytick_size_b = 0.5;                 % distance between y-ticks

%% Definitions phasor plot  -> modifiable
min_p = -2.0; max_p = 2.0;            % limits of phasor plot
xlabel_p = 'Re(F)';                   % xlabel in latex syntax
ylabel_p = 'Im(F)';                   % ylable in latex syntax
xtick_size_p = 0.5;                   % distance between x-ticks
ytick_size_p = 0.5;                   % distance between y-ticks
ind_c = [myBlue;myCyan;myYellow];     % color sequence of markers
mark = 'o';                           % shape of marker
mark_s = 9;                           % size of marker
mark_lw = 3;                          % linewidth of marker
quiv_sum = myRed;                     % color of resulting marker
headWidth = 5;                        % width of the individual phasor heads
headLength = 5;                       % length of the individual phasor heads

% Note, the individual phasors have always the same color as the according
% marker. The resulting phasor is twice als large and thick as the
% individual ones.

%% Defintions bode plot  -> modifiable
freq_c = myBlack;                     % color of bode plot
freq_lw = 1.5;                        % linewidth in bode plot
min_g = 0.0; max_g = 2.0;             % y-limits of gain

xlabel_f_g ='$\omega \tau_R$';        % xlabel gain in latex syntax 
xlabel_f_p = '$\omega \tau_R$';       % xlabel phase in latex syntax

ylabel_f_g = '$|F|$';                 % ylabel gain in latex syntax
ylabel_f_p = '$\angle F$';            % ylabel phase in latex syntax
xtick_size_f = 5;                     % distance between x-ticks
ytick_size_f_g = 0.5;                 % distance between y-ticks for gain
ytick_size_f_p = pi;                  % distance between y-ticks for phase

%% Definitions movie  -> modifiable
makeMovie = false;          % switch to "true" to create a movie. Warning, the code can take then a few minutes to run
movie_c = myBlue;           % color of marker and phasors in the movie
fps = 10;                   % framerate per second in movie. Warning, higher fps increases run time
stepsize = 8;               % frequency step size

%% ################ Everything below, DO NOT TOUCH! #######################
%%
%% ########################################################################
%% Computing the resulting transfer functions
tfs = repmat(h,1,length(omega)).*(exp(-1i*omega*tau'))'; % for every "h" the individual transfer functions are build
tf_sum = sum(tfs,1);                                     % sum of all individual transfer functions

%% Bar plot
f1=figure;
f1.Position = [50 50 1280 325];   % position and size of the whole figure (bars + phasors + bode) on the screen
subplot(2,3,[1,4])                % position of bar plot in the subplot arrangement

hold on
bar(tau/tau_R,h,bar_lw,'FaceColor',bar_c,'EdgeColor','none','LineWidth',1.5)      % plotting the bars
line([0 tau(end)/tau_R],[1 1],'LineStyle',':','LineWidth',0.5,'Color',[0.6 0.6 0.6]); % plotting dotted horizontal line at unity
set(gca,'box','off')              % no bounding box for bar plot
xlim([0 1])                       % x limits
ylim([min_b max_b])               % y limits  
xlabel(xlabel_b,'interpreter','latex','FontSize',16)  % x label
ylabel(ylabel_b,'interpreter','latex','FontSize',18)  % y label
xticks(0:xtick_size_b:1)                              % x ticks
yticks(min_b:ytick_size_b:max_b)                      % y ticks
hold off

%% Phasor plot
subplot(2,3,[2,5])
hold on
line([0 0],[min_p max_p],'LineStyle',':','LineWidth',0.5,'Color',[0.6 0.6 0.6]);                 % vertical dotted line
line([min_p max_p],[0 0],'LineStyle',':','LineWidth',0.5,'Color',[0.6 0.6 0.6]);                 % horizontal dotted line
plot(cos(0:0.001:2*pi),sin(0:0.001:2*pi),'LineStyle',':','LineWidth',0.5,'Color',[0.6 0.6 0.6]); % Unity circle
for i=1:1:length(ind) % plots all the markers
    plot(tf_sum(ind(i)),mark,'LineWidth',mark_lw,'Color',ind_c(i,:),'MarkerSize',mark_s);
end
% plots the phasor chaines for every marker
for j=1:1:length(ind) 
for i=1:1:length(h)
    % add only a phasor if "h" is non-zero
    if abs(h(i)) > eps  
    q = quiver(sum(real(tfs(1:i-1,ind(j))),1),sum(imag(tfs(1:i-1,ind(j))),1), ...
        real(tfs(i,ind(j))),imag(tfs(i,ind(j))),0,'MaxHeadSize',0.03,'Color',ind_c(j,:));
    % skip the first arrow for resizing the arrow heads
    if i>1 
    ah = annotation('arrow',...
            'headStyle','cback1','HeadLength',headLength,'HeadWidth',headWidth,'Color',ind_c(j,:));
        set(ah,'parent',gca);
        set(ah,'position',[q.XData q.YData 1*q.UData 1*q.VData]);
    end
    end
end
end
% plot the resulting phasor and resize its head
for j=1:1:length(ind) 
q_sum = quiver(0,0,real(tf_sum(ind(j))),imag(tf_sum(ind(j))),0,'MaxHeadSize',0.3,'LineWidth',2.5,'Color',quiv_sum);
ah = annotation('arrow',...
            'headStyle','cback1','HeadLength',2*headLength,'HeadWidth',2*headWidth,'Color',quiv_sum);
        set(ah,'parent',gca);
        set(ah,'position',[q_sum.XData q_sum.YData 1*q_sum.UData 1*q_sum.VData]);
end
%
% plot the transfer function in the complex plane
p = plot(tf_sum,'Color','k','LineWidth',1.5); 

xlabel(xlabel_p,'interpreter','latex','FontSize',16)
ylabel(ylabel_p,'interpreter','latex','FontSize',16)
ax3=gca;
ax3.XLim=[min_p max_p];
ax3.YLim=[min_p max_p];
xticks(min_p:xtick_size_p:max_p)
yticks(min_p:ytick_size_p:max_p)
ax3.Box='on'; % add bounding box for phasor plot
pbaspect([1 1 1]) % ensure that plot is quadratic and not strectched

%% Bode plot
% Gain
pos_g = [0.7 0.6 0.22 0.33];    % manually define position of the figure
subplot('Position',pos_g)       % position of bode GAIN in subplot arrangement
hold on
plot(omega*tau_R,abs(tf_sum),'LineWidth',freq_lw,'Color',freq_c)                  % plot gain curve
line([0 omega(end)],[1 1],'LineStyle',':','LineWidth',0.5,'Color',[0.6 0.6 0.6]); % dotted horizontal line
% plot all the markers
for i=1:1:length(ind)
b = plot(omega(ind(i))*tau_R,abs(tf_sum(ind(i))),mark,'LineWidth',mark_lw,'Color',ind_c(i,:),'MarkerSize',mark_s);
end
set(gca,'box','off')           % no bounding box for gain plot
xlim([0 omega(end)*tau_R])
ylim([min_g max_g])
xlabel(xlabel_f_g,'interpreter','latex','FontSize',16)
ylabel(ylabel_f_g,'interpreter','latex','FontSize',16)
hold off

% Phase
pos_p = [0.7 0.12 0.22 0.33];
subplot('Position',pos_p) % position of bode PHASE in subplot arrangement
hold on
plot(omega*tau_R,-1*wrapToPi(phase(tf_sum)),'LineWidth',freq_lw,'Color',freq_c)
% plot all the markers
for i=1:1:length(ind)
plot(omega(ind(i))*tau_R,-1*wrapToPi(phase(tf_sum(ind(i)))),mark,'LineWidth',mark_lw,'Color',ind_c(i,:),'MarkerSize',mark_s);
end
line([0 omega(end)],[pi pi],'LineStyle',':','LineWidth',0.5,'Color',[0.6 0.6 0.6]);   % dotted horizontal line at pi
line([0 omega(end)],[-pi -pi],'LineStyle',':','LineWidth',0.5,'Color',[0.6 0.6 0.6]); % dotted horizontal line at -pi
yticks([-pi 0 pi])
yticklabels({'-\pi','0','\pi'})
set(gca,'box','off') % no bounding box for phase plot
xlim([0 omega(end)]*tau_R)
xlabel(xlabel_f_p,'interpreter','latex','FontSize',16)
ylabel(ylabel_f_p,'interpreter','latex','FontSize',16)
hold off

%% Save figure as pdf
set(gcf, 'PaperPosition', [0 1 32 8]);  %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [32 9]);          %Set the paper to have width 5 and height 5.
saveas(gcf,name, 'pdf')
hold off


%% Make Movie
frames = size(length(tf_sum),2); % number of frames
ind = omega;
m = 1;    
if makeMovie

%The whole for loop below simply reproduces the above plot for all
%subsquent frequencies
for k=2:stepsize:length(omega)
f2=figure('Visible','off');
set(gcf,'color','w')
f2.Position = [50 50 1280 325]; % position and size of the whole figure (bars + phasors + bode) on the screen
subplot(2,3,[1,4])              % position of bar plot in the subplot arrangement
hold on
bar(tau/tau_R,h,bar_lw,'FaceColor',bar_c,'EdgeColor','none','LineWidth',1.5)      % plotting the bars
line([0 omega(end)],[1 1],'LineStyle',':','LineWidth',0.5,'Color',[0.6 0.6 0.6]); % plotting dotted horizontal line at unity
set(gca,'box','off') % no bounding box for bar plot
xlim([0 1]) 
ylim([min_b max_b]) 
xlabel(xlabel_b,'interpreter','latex')
ylabel(ylabel_b,'interpreter','latex')
xticks(0:xtick_size_b:1)
yticks(min_b:ytick_size_b:max_b)
hold off

%% Phasor plot
subplot(2,3,[2,5])
hold on
line([0 0],[min_p max_p],'LineStyle',':','LineWidth',0.5,'Color',[0.6 0.6 0.6]); % vertical dotted line
line([min_p max_p],[0 0],'LineStyle',':','LineWidth',0.5,'Color',[0.6 0.6 0.6]); % horizontal dotted line
plot(cos(0:0.001:2*pi),sin(0:0.001:2*pi),'LineStyle',':','LineWidth',0.5,'Color',[0.6 0.6 0.6]); % Einheitskreis
plot(tf_sum(k),mark,'LineWidth',mark_lw,'Color',movie_c,'MarkerSize',mark_s);
% plots the phasor chaines for every marker
for i=1:1:length(h)
    % add only a phasor if "h" is non-zero
    if abs(h(i)) > eps  
    q = quiver(sum(real(tfs(1:i-1,ind(k))),1),sum(imag(tfs(1:i-1,ind(k))),1), ...
        real(tfs(i,ind(k))),imag(tfs(i,ind(k))),0,'MaxHeadSize',0.03,'Color',movie_c);
    % skip the first arrow for resizing the arrow heads
    if i>1 
    ah = annotation('arrow',...
            'headStyle','cback1','HeadLength',headLength,'HeadWidth',headWidth,'Color',movie_c);
        set(ah,'parent',gca);
        set(ah,'position',[q.XData q.YData 1*q.UData 1*q.VData]);
    end
    end
end
% plot the resulting phasor and resize its head
q_sum = quiver(0,0,real(tf_sum(ind(k))),imag(tf_sum(ind(k))),0,'MaxHeadSize',0.3,'LineWidth',2.5,'Color',quiv_sum);
ah = annotation('arrow',...
            'headStyle','cback1','HeadLength',2*headLength,'HeadWidth',2*headWidth,'Color',quiv_sum);
        set(ah,'parent',gca);
        set(ah,'position',[q_sum.XData q_sum.YData 1*q_sum.UData 1*q_sum.VData]);
% plot the transfer function in the complex plane
p = plot(tf_sum,'Color','k','LineWidth',1.5); 

xlabel(xlabel_p,'interpreter','latex')
ylabel(ylabel_p,'interpreter','latex')
ax3=gca;
ax3.XLim=[min_p max_p];
ax3.YLim=[min_p max_p];
xticks(min_p:xtick_size_p:max_p)
yticks(min_p:ytick_size_p:max_p)
ax3.Box='on'; % add bounding box for phasor plot
pbaspect([1 1 1]) % ensure that plot is quadratic and not strectched

%% Bode plot
pos_g = [0.7 0.6 0.22 0.33];    % manually define position of the figure
subplot('Position',pos_g)       % position of bode GAIN in subplot arrangement
hold on
plot(omega*tau_R,abs(tf_sum),'LineWidth',freq_lw,'Color',freq_c) % plot gain curve
line([0 omega(end)],[1 1],'LineStyle',':','LineWidth',0.5,'Color',[0.6 0.6 0.6]); % dotted horizontal line
% plot all the markers
b = plot(omega(k)*tau_R,abs(tf_sum(k)),mark,'LineWidth',mark_lw,'Color',movie_c,'MarkerSize',mark_s);
set(gca,'box','off') % no bounding box for gain plot
xlim([0 omega(end)*tau_R])
ylim([min_g max_g])
xlabel(xlabel_f_g,'interpreter','latex')
ylabel(ylabel_f_g,'interpreter','latex')
hold off

pos_p = [0.7 0.12 0.22 0.33];
subplot('Position',pos_p) % position of bode PHASE in subplot arrangement
hold on
plot(omega*tau_R,-1*wrapToPi(phase(tf_sum)),'LineWidth',freq_lw,'Color',freq_c)
plot(omega(k)*tau_R,-1*wrapToPi(phase(tf_sum(k))),mark,'LineWidth',mark_lw,'Color',movie_c,'MarkerSize',mark_s);
line([0 omega(end)],[pi pi],'LineStyle',':','LineWidth',0.5,'Color',[0.6 0.6 0.6]); % dotted horizontal line at pi
line([0 omega(end)],[-pi -pi],'LineStyle',':','LineWidth',0.5,'Color',[0.6 0.6 0.6]); % dotted horizontal line at -pi
yticks([-pi 0 pi])
yticklabels({'-\pi','0','\pi'})
set(gca,'box','off') % no bounding box for phase plot
xlim([0 omega(end)]*tau_R)
xlabel(xlabel_f_p,'interpreter','latex')
ylabel(ylabel_f_p,'interpreter','latex')
hold off

M(m) = getframe(f2);

disp([int2str(m),' of ', int2str(ceil(length(omega)/stepsize)), ' frames'])
m = m+1;
    
end % for loop
%% 
      v = VideoWriter([pwd filesep name '_movie'], 'Uncompressed AVI');
      v.FrameRate=fps;
      open(v)
      writeVideo(v,M)
      close(v)
     
end % if makeMovie
