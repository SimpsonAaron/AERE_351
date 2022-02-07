clear; clc; close all
%%
%I) Given Parameters:
Alt_0 = 322; %km (Circular
R_e = 6378; %km
R_m = 1738; %km
R_s = 66183; %km R_SOI
D = 384400; %km
mu_e = 3.986*10^5; %km^3/s^2
mu_m = 4902.8; %km^3/s^2
omega_m = 2.6491*10^(-6); %rad/sec

%Converting to canonical:
D_DU = D/R_e; %DU
R_sDU = R_s/R_e; %DU
mu_eDU = 1;
%%
%%II) Test Departure Inputs
% %Project Test Inputs
% r_0 = 1.05; %DU
% v_0 = 1.372; %DU/TU
% phi_0 = 0; %deg
% lambda_1 = 30; %deg
% %[Alt_3,Em_0,Epsilon_2] = todamoon(r_0,v_0,phi_0,lambda_1);

% %Class Example 2 
% r_0 = 1.05; %DU
% v_0 = 1.38; %DU/TU
% phi_0 = 0; %deg
% lambda_1 = 26.5; %deg

%%
%III) Departure Conditions ()_0
r_0 = (Alt_0+R_e)/R_e; %DU
phi_0 = 0; %deg
v_0 =  sqrt(2/r_0)*.999; %DU/TU choosing 99.9% of escape velocity to get high speed and Em_0 < 0

%% Determining Lambda
lambda_1range = [25:0.000001:29]; %deg lambda range to test
Alt_3margin = .1; %
Alt_3plot = zeros(1,length(lambda_1range));
Epsilon_2plot = zeros(1,length(v_0ones));

for i = 1:length(lambda_1range)
    [Alt_3plot(i),Em_0,Epsilon_2plot(i)] = todamoon(r_0,v_0ones(i),phi_0,lambda_1range(i)); 
    %condensed function for cis-lunar trajectory at end of document
    if Em_0 < 0 & Epsilon_2plot(i) < 0 %"if approach is elliptical andretrograde"
        if abs(Alt_3plot(i)-90) < Alt_3margin %"if Periselenium altitude gets close to 90km, record valuse."
            lambda_1 = lambda_1range(i); 
            Alt_3margin = abs(Alt_3plot(i)-90);
            Alt_3 = Alt_3plot(i);
        end
    end
end

%We find that:
lambda_1; %deg
Alt_3;
Em_0;
plot (lambda_1range,Alt_3plot);
hold on ;
plot(lambda_1,Alt_3,'r*');
xlabel('Lambda_1');
ylabel('Altitude at Periselenium');
grid on;

%%
%IV) Arrival Outside of SOI:
Em_0 = v_0^2/2-mu_eDU/r_0;
Hm_0 = r_0*v_0*cosd(phi_0);
r_1 = sqrt(D_DU^2+R_sDU^2-2*D_DU*R_sDU*cosd(lambda_1)); %DU
v_1 = sqrt((Em_0+mu_eDU/r_1)*2); %equating energies
phi_1 = acosd(Hm_0/v_1/r_1);
gamma_1 = asind(R_sDU*sind(lambda_1)/r_1);
a_cl = mu_eDU/(2/r_0-v_0^2); %DU
e_cl = sqrt(1-Hm_0^2/a_cl/mu_eDU);
nu_1 = acosd(((a_cl*(1-e_cl^2))/r_1-1)/e_cl); %deg

%%
%V) Time of Flight:
E_1 = 2*atan(sqrt((1-e_cl)/(1+e_cl))*tan(deg2rad(nu_1)/2)) %rad
ToF_1 = sqrt(a_cl^3)*(E_1-e_cl*sin(E_1)); %TU;
ToF_1hrs = ToF_1*sqrt(R_e^3/mu_e)/360;0 %hrs

%%
%VI) Earth-Moon axis angle and departure angle:
alpha_1 = omega_m*ToF_1; %deg
gamma_0 = nu_1-alpha_1-gamma_1; %deg

%%
%VII) Inside SOI (SI units):
r_2 = R_sDU*R_e %km;
V_m = D*omega_m; %km/s
v_1= v_1*R_e/sqrt(R_e^3/mu_e); %km/s
v_2 = sqrt(V_m^2+v_1^2-2*V_m*v_1*cosd(phi_1-gamma_1)); %km/s
Epsilon_2 = asind(V_m/v_2*cosd(lambda_1)-v_1/v_2*cosd(lambda_1-(phi_1-gamma_1))); %deg
Em_SOI = v_2^2/2-mu_m/r_2 %km^2/s^2;
Hm_SOI = r_2*v_2*sind(Epsilon_2); %km/s
p_2 = Hm_SOI^2/mu_m; %km
a_2 = mu_m/(2*mu_m/r_2-v_2^2); %km
e_2 = sqrt(1-p_2/a_2);
r_3 = a_2*(1-e_2); %km

%%
%VIII) Time of Flight:
nu_2 = acosd((a_2*(1-e_2^2)/r_2-1)/e_2); %deg
F_2 = 2*atanh(sqrt((e_2-1)/(e_2+1))*tan(deg2rad(nu_2)/2)); %rad
ToF_2 = abs(sqrt(-a_2^3/mu_m)*(e_2*sinh(F_2)-F_2)); %s
ToF_2hrs = ToF_2/3600; %hrs
ToF_hrs = ToF_1hrs+ToF_2hrs; %hrs

%%
%IX) S/C Speed at Periselenium
v_3 = sqrt(2*(v_2^2/2-mu_m/r_2+mu_m/r_3)); %km/s
Alt_3 = r_3-R_m; %km

%%
%X) Answer Recap:
fprintf("v_0: %7.15f TU/DU",v_0);
fprintf("phi_0: %7.15f deg",phi_0);
fprintf("lambda_1: %7.15f deg",lambda_1);
fprintf("Em_0: %7.15f Du^2/TU^2",Em_0);
fprintf("nu_0: %7.15f deg",0);
fprintf("nu_1: %7.15f deg",nu_1);
fprintf("gamma_0: %7.15f deg",gamma_0);
fprintf("r_1: %7.15f DU",r_1);
fprintf("v_1: %7.15f TU/DU",v_1);
fprintf("phi_1: %7.15f deg",phi_1);
fprintf("gamma_1: %7.15f deg",gamma_1);
fprintf("a_cl: %7.15f DU",a_cl);
fprintf("e_cl: %7.15f",e_cl);
fprintf("v_2: %7.15f km/s",v_2);
fprintf("Epsilon_2: %7.15f deg",Epsilon_2);
fprintf("Em_2: %7.15f km^2/s^2",Em_SOI);
fprintf("a_2: %7.15f km",a_2);
fprintf("e_2: %7.15f",e_2);
fprintf("r_sp3: %7.15f km",r_3);
fprintf("Alt_3: %7.15f km",Alt_3);
fprintf("v_sp3: %7.15f km/s",v_3);


%%
%XI) Condensed Function for Plots:
function [Alt_3,Em_0,Epsilon_2] = todamoon(r_0,v_0,phi_0,lambda_1)
    Alt_0 = 322; %km (Circular
    R_e = 6378; %km
    R_m = 1738; %km
    R_s = 66183; %km
    D = 384400; %km
    mu_e = 3.986*10^5; %km^3/s^2
    mu_m = 4902.8; %km^3/s^2
    omega_m = 2.6491*10^(-6); %rad/sec
    
    %Converting to canonical:
    D_DU = D/R_e; %DU
    R_sDU = R_s/R_e; %DU
    mu_eDU = 1;
    
    Em_0 = v_0.^2./2-mu_eDU./r_0;
    Hm_0 = r_0*v_0*cosd(phi_0);
    r_1 = sqrt(D_DU^2+R_sDU^2-2*D_DU*R_sDU*cosd(lambda_1)); %DU
    v_1 = sqrt((Em_0+mu_eDU./r_1).*2); %equating energies
    phi_1 = acosd(Hm_0./v_1./r_1);
    gamma_1 = asind(R_sDU*sind(lambda_1)/r_1);
    
    a_cl = mu_eDU./(2./r_0-v_0.^2); %DU
    e_cl = sqrt(1-Hm_0.^2./a_cl./mu_eDU);
    nu_1 = acosd(((a_cl.*(1-e_cl.^2))./r_1-1)./e_cl); %deg
    
    E_1 = 2.*atan(sqrt((1-e_cl)./(1+e_cl)).*tan(deg2rad(nu_1)./2)); %rad
    ToF_1 = sqrt(a_cl.^3).*(E_1-e_cl.*sin(E_1)); %TU
    
    r_2 = R_sDU.*R_e; %km
    V_m = D.*omega_m; %km/s
    v_1= v_1.*R_e/sqrt(R_e^3./mu_e); %km/s
    v_2 = sqrt(V_m.^2+v_1.^2-2.*V_m.*v_1.*cosd(phi_1-gamma_1)); %km/s
    Epsilon_2 = asind(V_m./v_2.*cosd(lambda_1)-v_1./v_2.*cosd(lambda_1-(phi_1-gamma_1))); %deg
    Em_SOI = v_2.^2./2-mu_m./r_2; %km^2/s^2
    Hm_SOI = r_2.*v_2.*sind(Epsilon_2); %km/s
    p_2 = Hm_SOI.^2./mu_m; %km
    a_2 = mu_m./(2.*mu_m./r_2-v_2.^2); %km
    e_2 = sqrt(1-p_2./a_2);
    r_3 = a_2.*(1-e_2); %km
    
    nu_2 = acosd((a_2.*(1-e_2.^2)./r_2-1)./e_2); %deg
    F_2 = 2*atanh(sqrt((e_2-1)./(e_2+1)).*tan(deg2rad(nu_2)./2)); %rad
    ToF_2 = abs(sqrt(-a_2.^3./mu_m).*(e_2.*sinh(F_2)-F_2)); %s
    ToF = ToF_1*sqrt(R_e.^3./mu_e)+ToF_2; %s
    ToF_hrs = ToF./3600; %hrs
    
    v_3 = sqrt(2.*(v_2.^2./2-mu_m./r_2+mu_m./r_3)); %km/s
    Alt_3 = r_3-R_m; %km
end