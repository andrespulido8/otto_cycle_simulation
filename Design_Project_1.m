%% Design Project
%Author: Andres Pulido & Nicholas Mass
close all; clear; clc

%% Fuel constants
AF = 14.7; %Air-Fuel ratio []
Qhv = 43000; % Heating Value [kJ/kg]
Cp = 1.108; %Specific Heat at constant pressure [kJ/kg-K]
Cv = 0.821; %Specific Heat at constant volume [kJ/kg-K]
k = 1.35; %Average specific heat ratio []
R = 0.287; %Gas constant of air [kJ/kg-K]
nc = 1; %combustion efficiency

%% Given Constants
P1 = 100; %Inlet pressure [kPa]
rho = 1.255; %density of air at sea level [kg/m^3]
T1 = P1/(rho*R); %Ambient temperature [K]

%% 1st problem
rc = 1;
P3 = 0;
ntmax = 0;
while (P3 < 11000) %maximizing P3 will maximize nt for non-super
        rc = rc+0.01;
        T2 = T1*(rc)^(k-1);
        P2 = P1*(rc)^k;
        w12 = R*(T2 - T1)/(1-k);
        T3 = (Qhv*nc)/(AF+1)/Cv + T2;
        P3 = (P2/T2)*T3;
        T4 = T3*(1/rc)^(k-1);
        w34 = R*(T4-T3)/(1-k);
        wnet = w34 + w12;
        nt = 1-(1/rc)^(k-1); %Thermal Efficiency
        if (nt > ntmax) %record the largest net work and associated compression ratio
            ntmax = nt;
            rcnt = rc; %ideal rc to maximize imep for nonsupercharged engine
        end
end

fprintf("For problem 1:\n rc is %g and thermal efficiency is %g\n",rcnt,ntmax)

%% 2nd problem 
rc = 1;
P3 = 0;
imep = 0;
imepmax = 0;

while (P3 < 11000) 
        rc = rc+0.01;
        T2 = T1*(rc)^(k-1);
        P2 = P1*(rc)^k;
        w12 = R*(T2 - T1)/(1-k);
        T3 = ((Qhv*nc)/(AF+1)/Cv) + T2;
        P3 = (P2/T2)*T3;
        vtdc = Qhv/((AF+1)*Cv*(P3-P1*rc^k)); %calculate vtdc according to that rc
        T4 = T3*(1/rc)^(k-1);
        w34 = Cv*(T3-T4);
        wnet = w34 + w12;
        vbdc = rc*vtdc; %using previous vtdc but new rc
        imep = wnet/(vbdc - vtdc);
        if (imep > imepmax) %record the largest net work and associated compression ratio
            imepmax = imep;
            rcimep = rc; %ideal rc to maximize imep for nonsupercharged engine
        end
end

nt = 1-(1/rcimep)^(k-1); %Thermal Efficiency

fprintf("\nFor problem 2:\n imep is %g and thermal efficiency is %g for new rc equal to %g\n",imepmax,nt,rcimep)

%% figures 
v = linspace(vbdc,vtdc);
P = P1*(vbdc./v).^k; %first process

P = [P,linspace(P(100),P3)];
v(1,101:200) = v(100); %second process

P = [P,P(200)*(vtdc./flip(v(1:100))).^k]; 
v = [v,flip(v(1:100))]; %third process

P = [P,linspace(P(300),P(1),10)];
v(1,301:310) = v(300); %fourth process

figure;
plot(v,P,'LineWidth',1.5,'LineStyle','-','Marker','none','MarkerSize',15)
title('P vs v plot')
xlabel('specific volume [m^3/kg]')
ylabel('Pressure [kPa]')
ax = gca; %Get current axis
ax.FontName = 'Century Gothic'; ax.FontSize=16; ax.FontWeight='bold'; ax.GridAlpha=0.07; ax.GridLineStyle='--'; ax.LineWidth=1.5; ax.XColor=[0 0 0]; ax.XMinorTick='on'; ax.YColor=[0 0 0]; ax.YMinorTick='on';
