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

%% Constant given
P0 = 100; %Ambient pressure [kPa]
P1 = 100; %Inlet pressure to the engine with Supercharger [kPa]
rho = 1.255; %density of air at sea level [kg/m^3]
T0 = P1/(rho*R); %Ambient temperature [K]

%% 1st problem
P3 = 0;
P1nt = 0;
ntmax = 0;
rcnt = 0;
T5 = T0; 
P5 = P0; 
T6 = T0; 
P6 = P5; %p6 = p5 from looking at plot

while (P1 < 150)
    rc = 1;
    P1 = P1 + 1;
    P6a = P1; %p6a = p1 from looking at plot
    T1 = T0*(P1/P0)^((k-1)/k); %Output temperature of sc assuming isentropic efficiency is 100%
    while (P3 < 11000) 
            rc = rc+0.01;
            T2 = T1*(rc)^(k-1);
            P2 = P1*(rc)^k;
            w12 = R*(T2 - T1)/(1-k);
            T3 = (Qhv*nc)/((AF+1)*Cv) + T2;
            P3 = (P2/T2)*T3;
            vtdc = Qhv/((AF+1)*Cv*(P3-P1*rc^k)); %calculate vtdc according to that rc
            vbdc = rc*vtdc;
            T4 = T3*(1/rc)^(k-1);
            w34 = R*(T4-T3)/(1-k);
            wsc = Cp*(T1-T0); 
            %lower loop
            v6a = (vbdc-vtdc)/3;
            wlowerloop = (vbdc - v6a)*(P1 - P5) + 0.5*(v6a - vtdc)*(P1 - P6);
            wnet = w12 + w34 + wlowerloop - wsc;
            qin = Cv*(T3-T2);
            nt = wnet/qin;
            if (nt > ntmax) %record the largest net work and associated compression ratio
                ntmax = nt;
                rcnt = rc; %ideal rc to maximize nt for nonsupercharged engine
                P1nt = P1;
            end
    end
end

fprintf("For problem 1:\n rc is %g at an inlet pressure of %g and thermal efficiency is %g\n",rcnt,P1nt,nt)

%% 2nd problem 
rc = 1;
P3 = 0;
P1 = 100;
imepmax = 0;
rcimep = 0;
wnet = 0;

while (P1 < 150)
    P1 = P1 + 1;
    P6a = P1; %p6a = p1 from looking at plot
    T1 = T0*(P1/P0)^((k-1)/k); %Output temperature of sc assuming isentropic efficiency is 100%
    while (P3 < 11000) 
        rc = rc+0.01;
        T2 = T1*(rc)^(k-1);
        P2 = P1*(rc)^k;
        w12 = R*(T2 - T1)/(1-k);
        T3 = ((Qhv*nc)/(AF+1)/Cv) + T2;
        P3 = (P2/T2)*T3;
        vtdc = Qhv/((AF+1)*Cv*(P3-P1*rc^k)); %calculate vtdc according to that rc
        vbdc = rc*vtdc;
        T4 = T3*(1/rc)^(k-1);
        w34 = R*(T4-T3)/(1-k);
        wsc = Cp*(T1-T0); %work needed to drive the supercharger
        % lower loop
        v6a = (vbdc-vtdc)/3;
        wlowerloop = (vbdc - v6a)*(P1 - P5) + 0.5*(v6a - vtdc)*(P1 - P6);
        wnet = w12 + w34 + wlowerloop - wsc;
        vbdc = rc*vtdc; %using previous vtdc but new rc
        imep = wnet/(vbdc - vtdc);
        if (imep > imepmax) %record the largest net work and associated compression ratio
            imepmax = imep;
            rcimep = rc; %ideal rc to maximize imep for nonsupercharged engine
            wnetmax = wnet;
            qin = Cv*(T3-T2);
            P1imep = P1;
        end
    end
end
ntimep = wnetmax/qin;

fprintf("\nFor problem 2:\n imep is %g and thermal efficiency is %g for new rc equal to %g at inlet pressure of %g\n",imepmax,ntimep,rcimep,P1imep)

%% figures 
v = linspace(vbdc,vtdc);
P = P1*(vbdc./v).^k; %first process

P = [P,linspace(P(100),P3)];
v(1,101:200) = v(100); %second process

P = [P,P(200)*(vtdc./flip(v(1:100))).^k]; 
v = [v,flip(v(1:100))]; %third process

P = [P,linspace(P(300),P0)];
v(1,301:400) = v(300); %fourth process

P(1,401:500) = P(400);
v = [v,linspace(vbdc,vtdc)];

P = [P,linspace(P0,P1,50)];
v = [v,linspace(vtdc,v6a,50)];

P(1,551:600) = P(550);
v = [v,linspace(v6a,vbdc,50)];

figure;
plot(v,P,'LineWidth',1.5,'LineStyle','-','Marker','none','MarkerSize',15)
title('P vs v plot')
xlabel('specific volume [m^3/kg]')
ylabel('Pressure [kPa]')
ax = gca; %Get current axis
ax.FontName = 'Century Gothic'; ax.FontSize=16; ax.FontWeight='bold'; ax.GridAlpha=0.07; ax.GridLineStyle='--'; ax.LineWidth=1.5; ax.XColor=[0 0 0]; ax.XMinorTick='on'; ax.YColor=[0 0 0]; ax.YMinorTick='on';
