%% Function to compute the X and Y matrices:
close all
clear all
clc

syms M1 M2 rho1 rho2 c1 c2 A1 A2 p1 p2 T1 T2 u1 u2 Cp
syms ft ft_Tu gt gt_Tu jt jt_Td ht ht_Td st st_Td Qh V
% Inlet
ft = (1-M1)/(1+M1)*gt_Tu + V*c1/(A1*(1 + M1));

% ft = -(1-M1)^2/(1+M1)^2*gt_Tu + V/(1 + M1)^2;


% outlet
jt = -ht_Td;
% Jumps
m1 = A1*rho1*u1;
m2 = A2*rho2*u2;
H1 = Cp*T1 + 1/2*u1^2;
H2 = Cp*T2 + 1/2*u2^2;

ph1 = ft + gt;
ph2 = ht + jt;
rh1 = 1/(c1^2)*(ft + gt);
% rh2 = 1/c2^2*(ht + jt + st);
uh1 = 1/(rho1*c1)*(ft - gt);
uh2 = 1/(rho2*c2)*(ht - jt);

mh1 = A1*(rho1*uh1 + rh1*u1);


rh2 = rho2*(mh1/m1 - uh2/u2);


mh2 = A2*(rho2*uh2 + rh2*u2);





fh1 = A1*ph1 + mh1*u1 + m1*uh1;
fh2 = A2*ph2 + mh2*u2 + m2*uh2; % Subs mh2 with mh1

Th1 = T1*(ph1/p1 - rh1/rho1);
Th2 = T2*(ph2/p2 - rh2/rho2);

Hh1 = Cp*Th1 + u1*uh1;
Hh2 = Cp*Th2 + u2*uh2;

eh1 = mh1*H1 + m1*Hh1;
eh2 = mh2*H2 + m2*Hh2; % Subs mh2 with mh1

% syms gamma R_gas
% eh1x = A1*( (gamma/(gamma-1)*(ph1*u1 + p1*uh1) + rh1*u1*1/2*u1^2 + rho1*uh1*1/2*u1^2 + rho1*u1*u1*uh1));
% eh2x = A2*( (gamma/(gamma-1)*(ph2*u2 + p2*uh2) + rh1*u1*1/2*u2^2 + rho1*uh1*1/2*u2^2 + rho1*u1*u2*uh2));

% 
% % Whay are they different:
% SS = (eh2 - eh2x);
% 
% SS = subs(SS,[T1,T2,Cp],[p1/(rho1*R_gas),p2/(rho2*R_gas),gamma*R_gas/(gamma-1)]);
% SS = subs(SS,[M1,M2],[u1/c1,u2/c2]);
% SS = subs(SS,[p1,p2],[rho1*c1^2/gamma,rho2*c2^2/gamma]);
% SS = subs(SS,u2,m1/(rho2*A2));
% % SS = subs(SS,A1,A2)
% SS = simplify(SS)


% return


E4 = fh2 - fh1 - ((A2 - A1)*ph1);
E5 = eh2 - eh1 - Qh;
% Set equations
Eq = [E4 E5/c1].';

% Aim is to get expressions of the form X*ft - Y*ft_T - F = 0

X(:,1) = diff(Eq,gt);
X(:,2) = diff(Eq,ht);

Y(:,1) = diff(Eq,gt_Tu);
Y(:,2) = diff(Eq,ht_Td);

Y = -Y;

F = -Eq + (X*[gt;ht] - Y*[gt_Tu;ht_Td]);
F = simplify(F);

% Simplify
W = [X Y F];
syms gamma R_gas m
W = W/A1;
W = subs(W,[T1,T2,Cp],[p1/(rho1*R_gas),p2/(rho2*R_gas),gamma*R_gas/(gamma-1)]);
W = subs(W,[p1,p2],[rho1*c1^2/gamma,rho2*c2^2/gamma]);
W = subs(W,[rho2,rho1],[m/(A2*u2),m/(A1*u1)]);
W = subs(W,[u1,u2],[M1*c1,M2*c2]);

W = simplify(W,1000);

Xr = W(:,1:2);
Yr = W(:,3:4);
F  = W(:,end);

XX = Xr;
for j = 1:2
    for k = 1:2
       fprintf(['X(%d,%d) = ',char(XX(j,k)),'; \n'],[j,k])
    end
end

XX = Yr;
for j = 1:2
    for k = 1:2
       fprintf(['Y(%d,%d) = ',char(XX(j,k)),'; \n'],[j,k])
    end
end



% SS = eh1/A1 - (gamma/(gamma-1)*(ph1*u1 + p1*uh1) + rh1*1/2*u1^3 + rho1*3/2*u1^2*uh1);
% % SS = eh2/A2 - (gamma/(gamma-1)*(ph2*u2 + p2*uh2) + rh2*1/2*u2^3 + rho2*3/2*u2^2*uh2);
% SS = subs(SS,[T1,T2,Cp],[p1/(rho1*R_gas),p2/(rho2*R_gas),gamma*R_gas/(gamma-1)]);
% SS = subs(SS,[u1,u2],[M1*c1,M2*c2]);
% SS = subs(SS,[p1,p2],[rho1*c1^2/gamma,rho2*c2^2/gamma]);
% SS = simplify(SS)
