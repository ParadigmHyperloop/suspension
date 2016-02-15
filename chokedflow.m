%% Worst case choked flow around skates 
% Prepared for OpenLoop Hyperloop design team
% for the SpaceX hyperloop design competition
% Updated 02/14/16
clear;
 
% the purpose of this script is to calculate the worst case mass flow rate
% to generate a specified ride height given choked flow at air skate
% boundaries

% Air temperature based on adiabatic expansion of tank air at 300 Kelvin
%T = 261;       % [K] from 3000 psi tanks
T = 275.3;       % [K] from 1560 psi tanks
R = 8.3144598;  % [J/mol*K] molar gas constant
MW = 0.02897;   % [kg/mol] molecular weight of air
Rs = R/MW;      % [J/kg*K] specific gas constant of air
k = 1.40;       % ratio of specific heats
 
p_amb = 1000;         % [Pa] ambient pressure
a = sqrt(k*Rs*T);     % speed of sound [m/s]

% calculate area of boundary
height = .003;          % [m]
A = 15.8496 * height;   % [m^2]  2* perimeter * height (geometry of skates)

% calculate skate area
area_skates_Total = 2.22967; % [m^2]
 
% pod mass
mass = 750;    % [kg]

g = 9.81;        % [m/s^2]
p = mass / area_skates_Total * g;  % gauge pressure at skates to support pod mass
p = p+ p_amb;  % absolute pressure at skates
rho = p/Rs/T;   % density of air
 
mdot = rho* A * a; % [kg/s]  mass flow if gas with density rho flows through area A

fprintf('\nFor a ride gap height of %.1f mm\n' ,height*1e3);
fprintf('Maximum mass flow rate=  %.3f kg/s\n', mdot);
fprintf('At a pressure of %.0f Pa\n', p);
fprintf('At an ambient pressure of %.0f Pa\n', p_amb);