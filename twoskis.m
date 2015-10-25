% Simulates a pair of long, rectangular airskates

% Skate parameters
k = 8e-8*(2.54e-2)^2;		% Air permeability [m^2]
D = 0.187*2.54e-2;			% Thickness of porous layer [m]
P0 = 11e3;					% Internal pressure of the skate [Pa]
n = 2;						% Number of skates
W = 0.3048;					% Skate width [m]
L = 10*0.3048;				% Skate length [m]
T = 400;					% Nominal temperature of air through skate [K]

% Independent variable - Gap height
H = linspace(0,1500e-6,1e3); % Gap height [m]

% Other system parameters
P_supp = 2.1e3;				% Pressure of supply air to skate compressor [Pa]
T_supp = 300;				% Temperature of supply air to skate compressor [K]
P_tube = 99;				% Tube pressure [Pa]
m_pod = [500 1000 1500];	% Pod mass [kg]

% Physical constants
M_air = 28.97e-3;			% Molecular weight of air [kg/mol]
R = 8.3144598;				% Molar gas constant [J/K*mol]
g = 9.81;					% Acceleration of gravity [m/s^2]
c = 343;					% Speed of sound [m/s]

% Intermediate calculations
mu = 0.01827e-3*(291.15+120)/(T+120)... % Viscosity of air [Pa*s]
	*(T/291.15)^(1.5);
alpha =	sqrt(12*k./(H.^3*D));	% Dimensionless parameter "alpha"
A = n*L*W;						% Total skate area [m^2]

% Force as a function of gap height [N]
F = n*L*P0*(W-2./alpha.*tanh(alpha.*W/2));

% Flow rate as a function of gap height [kg/s]
m_flow = (P0*M_air/(R*T))*...
	W*P0*alpha./(2*mu).*tanh(alpha*W/2).*H.^3.*(1/2-1/3).*(2*(W+L));

% Derivative of force as a function of gap height [N/m]
F_H = n*L*P0*(2./alpha.^2.*tanh(alpha.*W/2) ...
			  - W./alpha.*sech(alpha.*W/2).^2) ...
	  .*alpha.*(-3/2).*(1./sqrt(H));

% Calculate frequency of oscillation
m_pod_H = F/g;		% Mass as a function of equilibrium ride height [kg]
freq_H = 1/(2*pi).*sqrt(-(1./m_pod_H).*F_H); % Frequency of oscillation as
											 % a function of height [Hz]

% Plot carrying force vs gap height
figure(1)
subplot(2,1,1)
plot(H*1e6,F/g)
ylabel('Force [kg_F]')
xlabel('Gap Height [\mum]')

% Plot airflow as a function of gap height
subplot(2,1,2)
plot(H*1e6,m_flow)
ylabel('Flow Rate [kg/s]')
xlabel('Gap Height [\mum]')

% Plot vibration frequency as a function of equilibrium gap height
figure(2)
subplot(2,1,1)
plot(H*1e6,freq_H)
ylabel('Vibration Freq. [Hz]')
xlabel('Equilibrium Gap Height [\mum]')

% Plot vibration frequency as a function of pod mass
subplot(2,1,2)
plot(m_pod_H,freq_H)
ylabel('Vibration Freq. [Hz]')
xlabel('Pod Mass [kg]')

% Determine temperature of the uncooled air based on isentropic compression
V1 = 1;
V2 = (P_supp/P0)^(1/1.4);
T_uncooled = P0*V2/(P_supp*V1/T);	% Uncooled air temperature [K]

fprintf('%i skates at %g m^2 per skate\n',n,pi*R^2)
fprintf('Total skate area: %g m^2\n',A)
fprintf('Uncooled air temperature: %g K\n',T_uncooled)