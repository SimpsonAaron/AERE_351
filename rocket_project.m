clear, clc, close all

mu = 3.986*10^5;

nu = degtorad(25); %25 degrees in rad
rP = 72+6378; %r perigee GTO
rA = 35786+6378; %r apogee GTO

a = (rA+rP)/2; %semi-major axis (km) GTO
e = (rA-rP)/(rA+rP); %eccentricity GTO
p = a*(1-e^2); 

r = p/(1+e*cos(nu)); %r at nu (km) GTO
v = sqrt(2*mu/r-mu/a); %velocity at nu (km/s) GTO

v_surface = 0.4084; %earth surface velocity (km/s)
v_lost = 1.5; %extra velocity lost to drag, gravity, etc. (km/s)

V_GTO = v + v_lost - v_surface; %Final velocity for GTO insertion at nu (km/s)


V_GTO_a = sqrt(2*mu/rA-mu/a); %velocity at apogee of transfer orbit (km/s)
V_c2 = sqrt(mu/(rA)); %velocity of circular geostationary orbit 
i=deg2rad(28.5); %change in inclination
V3 = sqrt(V_GTO_a^2+V_c2^2-2*V_c2*V_GTO_a*cos(i)); %delta v (burnout) needed to get 3rd stage to Geostationary from GTO
g = 9.81*10^-3; %km/s^2 earth gravity constant
m_03 = 3200; %total stage 3 weight (kg)
Isp3 = 240; %(sec) stage 3 specific impulse
E3 = 0.15; %stage 3 structure ratio
V_ex3 = Isp3*g; %stage 3 exhaust velocity (km/s)
pi_3 = (exp(V3/-V_ex3)-E3)/(1-E3); %stage 3 structure ratio
m_star3 = pi_3*m_03; %stage 3 electronic payload (kg)


pi_2_range = linspace(.001,1);
m_0_output = totalmass(pi_2_range,V_GTO); %using the function defined at the end of the script
plot(pi_2_range,m_0_output,'r');
title('Total Mass as a Function of 2nd Stage Payload Ratio');
xlabel('pi_2');
ylabel('m_0');
grid on


x = .0001:0.00001:.2;
y = totalmass(x,V_GTO); %using the function defined below
idx = islocalmin(y);
figure(1)
hold on
plot(x(idx),y(idx),'*b')
legend('Curve','Local Min')

hold off
fprintf('Min located at %0.5f\n',x(idx))


Isp1 = 280; %sec 1st stage specific impulse
Isp2 = 455; %sec 2nd stage specific impulse

E1 = 0.11; %1st stage structure ratio
E2 = 0.13; %2nd stage structure ratio

V_ex1 = Isp1*g; %1st stage exhaust velocity
V_ex2 = Isp2*g; %2nd stage exhaust velocity



pi_2 = x(idx); %input payload ratio value based on minimum calculated above
m_02 = m_03./pi_2; %stage 2 total mass (kg)
m_s2 = (m_03./pi_2-m_03)*E2; %stage 2 structure mass (kg)
m_p2 = m_02-m_03-m_s2; %stage 2 propellant mass (kg)
V2 = -V_ex2*log(E2+(1-E2)*pi_2); %stage 2 burnout velocity


V1 = V_GTO-V2; %stage 1 burnout velocity
m_01 = totalmass(pi_2,V_GTO);% total/stage1 1 mass (kg)
m_s1 = (m_01-m_02)*E1; %stage 1 structure mass (kg)
m_p1 = m_01-m_02-m_s1; %stage 1 propellant mass (kg)
pi_1 = m_02/m_01; %stage 1 payload ratio 

function [m_0] = totalmass(pi_2,V_GTO); %this function outputs m_0 as a function of the pi_2 and V_GTO inputs
g = 9.81*10^-3; %km/s^2 earth gravity constant

Isp1 = 280; %sec 1st stage specific impulse
Isp2 = 455; %sec 2nd stage specific impulse

E1 = 0.11; %1st stage structure ratio
E2 = 0.13; %2nd stage structure ratio

V_ex1 = Isp1*g; %1st stage exhaust velocity
V_ex2 = Isp2*g; %2nd stage exhaust velocity

m_03 = 3200; %3rd stage mass (kg)

%%stage 2 masses based on pi_2:
m_02 = m_03./pi_2;
m_s2 = (m_03./pi_2-m_03)*E2;
m_p2 = m_02-m_03-m_s2;

%%stage 1 masses based on pi_2
m_0 = ((1-E1).*(m_03./pi_2))./(((exp(V_GTO./-V_ex1))./((E2+(1-E2).*pi_2).^(V_ex2./V_ex1)))-E1); %Total rocket mass as a function of pi2
end