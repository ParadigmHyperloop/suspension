% Testing the effect of different tube pressures
% on a pair of long, rectangular airskates

% Physical constants
M_air = 28.97e-3;     % Molecular weight of air [kg/mol]
R = 8.3144598;        % Molar gas constant [J/K*mol]
g = 9.81;             % Acceleration of gravity [m/s^2]
c = 343;              % Speed of sound [m/s]

% Skate parameters
k = 5e-11;              % Air permeability [m^2]
D = 5e-3;               % Thickness of porous layer [m]
n = 2;                  % Number of skates
W = 0.3048;             % Skate width [m]
L = 10*0.3048;          % Skate length [m]
T = 400;                % Nominal temperature of air through skate [K]
H = 0.75*1e-3;			% Nominal gap height [m]
% Other system parameters
P_supp = 2.1e3;        % Pressure of supply air to skate compressor [Pa]
T_supp = 300;          % Temperature of supply air to skate compressor [K]
m_pod = 3000;  % Pod mass [kg]

% Independent variable - Tube Pressure
P_tube = 6894.76*(0.02:0.01:4.0); % Tube pressure [Pa]

% Internal skate air pressure set to provide twice the weight of the pod
P0 = m_pod*g/(n*W*L) + P_tube; % Internal skate pressure [Pa]

% Intermediate calculations
P0_gauge = P0 - P_tube;
mu = 0.01827e-3*(291.15+120)/(T+120)... % Viscosity of air [Pa*s]
  *(T/291.15)^(1.5);
alpha = sqrt(12*k./(H.^3*D));          % Dimensionless parameter "alpha"
A = n*L*W;                              % Total skate area [m^2]

% Required skate pressure as a function of tube pressure [Pa]
P0_needed = (m_pod*g)./(n*L*(W-2./alpha.*tanh(alpha.*W/2))) + P_tube;

% Flow rate as a function of tube pressure [Pa]
m_flow = (M_air/(R*T)).*...
  W*P0.*P0_gauge.*alpha./(2*mu).*tanh(alpha*W/2).*H.^3.*(1/2-1/3).*(2*(W+L));

% Derivative of force as a function of gap height [N/m]
F_H = n*L*P0_gauge.*(2./alpha.^2.*tanh(alpha.*W/2) ...
        - W./alpha.*sech(alpha.*W/2).^2) ...
    .*alpha.*(-3/2).*(1./sqrt(H));

% Calculate frequency of oscillation
freq = 1/(2*pi).*sqrt(-(1./m_pod).*F_H);  % Frequency of oscillation as a
                                          % function of tube pressure [Hz]

% Plot required skate pressure vs tube pressure
figure(1)
subplot(3,1,1)
plot(P_tube/100,P0_needed/100)
ylabel('Skate pressure [mbar]')
xlabel('Tube Pressure [mbar]')

% Plot airflow as a function of tube pressure
subplot(3,1,2)
plot(P_tube/100,m_flow)
ylabel('Flow Rate [kg/s]')
xlabel('Tube Pressure [mbar]')

% Plot frequency of oscillation as a function of tube pressure
subplot(3,1,3)
plot(P_tube/100,freq)
ylabel('Frequency [Hz]')
xlabel('Tube Pressure [mbar]')
