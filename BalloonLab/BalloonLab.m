close all; clear; clc;

%% variable declaration
m_payload = 500; %kg 

rho_hydrogen_day = 0.00226; %kg/m^3
rho_hydrogen_night = 0.00343; %kg/m^3

rho_air = .04008; %kg/m^3 -- at 25km altitude   
rho_hydrogen = 0.0027; %kg/m^3 -- at 25km altitude -- hydrogen
rho_mat = 920; %kg/m^3

FS = 1.5; %unitless
YS = 123.44 * 1e6; %Pa

%% calculating specific values 
air = ((4*pi) / 3) * rho_air;
gas = ((4*pi) / 3) * rho_hydrogen;
mat = (((40*pi) * (1.5)) / (2 * YS)) * rho_mat; 

%% calculate value for radius
radius = nthroot((m_payload / (air - gas - mat)), 3); %m 

%% Calculate volume (pre thermal effects)
V = (4/3) * pi * (radius)^3; %m^3

%% Calculate volume (thermal effects included)
[volume_day, volume_night, radius_day, radius_night] = volumes(rho_hydrogen, V, rho_hydrogen_day, rho_hydrogen_night);

%% Thickness
thickness = ((FS) * 10 * radius) / (2 * YS) * 1e6; %mu m

%% plots
Factor_safety = 0:.5:25;
thick = ((Factor_safety * 10 * radius) / (2 * YS)) * 1e6; %mu m
p = plot(Factor_safety, thick);
xlabel('Factor of Safety')
ylabel('Balloon Thickness [$\mu$m]')
title('Factor of Safety vs Balloon Thickness')
hold on 
q = plot(FS, thickness, 'ro');
legend('Thickness', 'Selected Thickness and FOS')

set(0,'defaultTextInterpreter','latex')
set(gca,'FontSize',18)
p(1).LineWidth = 3;
p(1).MarkerSize = 15;
q(1).LineWidth = 5;
xlim([1 10])
grid on

%% Air density
%day
a_d = (40 * pi * radius_day^3 * FS) / (2 * YS);
b_d = (4/3)*(pi * radius_day^3 * rho_hydrogen_day);
c_d = 4*pi*radius_day^3;

rho_air_day = 3 * (600 + a_d + b_d) / (c_d);  %balast = 100kg

%night
a_n = (40 * pi * radius_night^3 * FS) / (2 * YS);
b_n = (4/3)*(pi * radius_night^3 * rho_hydrogen_night);
c_n = 4*pi*radius_night^3;

rho_air_night = 3 * (500 + a_n + b_n) / (c_n); 

%% functions
function [volume_day, volume_night, radius_day, radius_night] = volumes(rho_gas, V0, rho_gas_day, rho_gas_night)

    mass_H = 40.5 %rho_gas * V0 orignially equaled 36.5 kg - altitude was too low 

%% day
    volume_day = mass_H / rho_gas_day; %m^3
    radius_day = nthroot(((3 * volume_day) / (4*pi)), 3);
    
 %% Night
    volume_night = mass_H / rho_gas_night; %m^3
    radius_night = nthroot(((3 * volume_night) / (4*pi)), 3);
    
end