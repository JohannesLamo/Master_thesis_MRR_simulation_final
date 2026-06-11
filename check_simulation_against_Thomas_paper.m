clear all; 
clc; 

% This script compares my simulation results to figures from the papers: 
% - "Cascadable excitability in microrings", by "Thomas Van Vaerenbergh"
% - "Simplified description of self-pulsation and excitability by thermal
% and free-carrier effects in semiconductor microcavities ", by "Thomas Van Vaerenbergh"



%% 1) Generation of data sequence

Ltrain = 800;                                            % Length of Training set
Ltest = 0;                                              % Length of Testing set
Warm1 = 0;                                              % Length of NARMA-10 sequence warm-up
Warm2 = 0;                                              % Length of the reservoir warm-up. 
                                                        % We also use this sequence between training and test
                                                        % To avoid memory of
                                                        % training data in
                                                        % MRR when running
                                                        % test.

tot_Ltrain = Ltrain + Warm2;                            % Total number of bits for training = Ltrain + warm up
tot_Ltest = Ltest + Warm2;                              % Total number of bits for test = Ltest + warm up
total_L   = tot_Ltrain + tot_Ltest;                     % Total number of bits for both training and testing. 



% Create mask and set number of time steps per bit. 
% Since we work with a discrete time grid, it is important that both N_mask_raw, steps_per_bit_raw and
% (steps_per_bit_raw / N_mask_raw) are integers. If (steps_per_bit_raw /
% N_mask_raw) is not an integer, the function "change_N_mask" will alter the
% values of
% either N_mask_raw, steps_per_bit or both in order to fulfill the demand.
N_mask_raw = 50;                                         % Number of wanted virtal nodes.
steps_per_bit_raw = 500;                                 % Time steps per bit. Must be an integer.
[N_mask, steps_per_bit] = change_N_mask(steps_per_bit_raw, N_mask_raw); 

% Print message if N_mask or steps_per_bit have been altered. 
if N_mask == 1
   disp(['steps_per_bit / N_mask needs to be an integer. Since this is not the case, N_mask is set to N_mask = 1. The data is therefore not masked.' ...
       ' Change steps_per_bit.'])
elseif N_mask ~= N_mask_raw || steps_per_bit ~= steps_per_bit_raw
    disp(['steps_per_bit / N_mask_raw is not an integer. Therefore, N_mask is changed to ', num2str(N_mask) ', and steps_per_bit is changed by ', num2str(steps_per_bit - steps_per_bit_raw) ' number of steps per bit.'])
end    


bitlength = 1e-9;                                       % Bit width [s]
solver_steps = bitlength/steps_per_bit;                 % Solver steps in RK4 methods [s]
time_window = total_L*bitlength;                        % Total time window for our data sequence [s]
time_axis = 0:solver_steps:(time_window-solver_steps);  % Time axis [s]
bias = 8;


%% 2) Parameters for Micro Ring Resonator. Same as in paper "Cascadable excitability in microrings"

nsi = 3.476;                                            % Unpertubated refractive index for Silicon [~]
cp = 0.7e3;                                             % specific heat capacity at constant pressure [J*(kg*K)^(-1)]

rho = 2330;
Vth = 3.19e-18;
m = rho*Vth;

L = 2*pi*4.0e-6;                                        % Circumerence MRR [m]
B_TPA = 8.4e-12;                                        % TPA coefficient [m*W^(-1)]
alpha_dB = 0.8*1e2;                                     % Attenuation per meter in desibel [dB/m]
alpha_lin = alpha_dB*log(10)/10;                        % Attenuation per meter linear [1/m]
zeta = 0.95;                                            % Coupling coefficient over to delay loop
phi_d = 0;                                              % Phase shift at Add-port.

dnsi_dT = 1.86e-4;                                      % d(nsi)/dT [K^(-1)]
dnsi_dN = -1.73e-27;                                    % d(nsi)/dN [m^(-3)]

GAMMA_TH = 0.9355;                                      % Confinement factor for temperature [~]
GAMMA_TPA = 0.9964;
GAMMA_FCA = 0.9996;                                     % FCA confinement factor [~]
sigma_fca = 1e-21;                                      % FCA coefficient in Silicon [m^2]
V_FCA = 2.36e-18;                                       % FCA effective mode volume [m^3]
V_TPA = 2.59e-18;                                       % TPA effective mode volume [m^3]

c = 299792458;                                          % Speed of light [m/s]
h = 6.62607015e-34 ;                                    % Plancks constant [J*s]
h_bar = h/2/pi;                                         % Reduced Plancks constant [J*s]
lambda_resonance_0 = 1552.770e-9;                        % Resonance wavelength in cold cavity [m]
omega_resonance_0 = 2*pi*c/lambda_resonance_0;          % Resonance angular frequency in cold cavity [rad/s]


% Time constants
tau_c = 205e-12 ;                                       % Coupling lifetime [s]. ALL-PASS CONFIG
tau_r = tau_c; 
tau_ph = (1/tau_r + 1/tau_c)^-1;                        % Total loss lifetime. Should be shorter than the shortest of the two parameters

%%%% DESSE ER ENDRA PÅ
tau_th = 65e-9; 
tau_fc = 5.3e-9; 

tau_d = 0.5e-9; 
t_sq = 0.9604;                                          % Self-coupling power coeffifient (t^2)
k_sq = 1-t_sq;                                          % Cross-couplung power coefficient. (k^2). Remember: t_sq + k_sq = 1
eta_lin = 0.4;                                            % Fraction of linear loss due to absorbtion. 



%% Test RK4 solver
mode = "all_pass";

%% USE this if comparing to fig 3. in the paper "Cascadable excitability in microrings"
delta_lambda = 62e-12; 
omega_p = 2*pi*c / (delta_lambda + 2*pi*c/omega_resonance_0);
Pin = 0.6e-3;
Ein = sqrt(Pin)*ones(size(time_axis));
% Initial conditions. SI-units
T_0 = 0.7; 
N_0 = 0e22; 
a_0 = 0;

%% USE THIS if comparing to figure 1 in the paper "Simplified description of self-pulsation and excitability by thermal and free-carrier effects in semiconductor microcavities "
delta = -3; % From paper 
omega_p = delta/tau_ph + omega_resonance_0; 
Pin_avg_lin = 0.6e-3;
Ein = sqrt(Pin_avg_lin)*ones(size(time_axis));
% Initial conditions. SI-units
T_0 = 0.0; 
N_0 = 0e22; 
a_0 = 0;


%% 4D version. USE THIS if comparing to the paper "Cascadable excitability in microrings"

[E_drop, E_through, a_norm, theta_norm, n_norm, tau, a, delta_Temp, delta_N] = RK4_4D_norm(time_axis, Ein, T_0, N_0, a_0, solver_steps, tau_d, ...
    tau_th, tau_fc, tau_r, tau_c, eta_lin, omega_resonance_0, omega_p, m, cp, nsi, ...
    dnsi_dT, dnsi_dN, GAMMA_TH, GAMMA_FCA, GAMMA_TPA, V_FCA, V_TPA, B_TPA, c, h_bar, sigma_fca, zeta, phi_d, mode);


%% 2D version. USE THIS if comparing to the paper "implified description of self-pulsation and excitability by thermal and free-carrier effects in semiconductor microcavities"

[E_drop_2D, E_through_2D, a_norm_2D, theta_norm_2D, n_norm_2D, tau_2D, a_2D, delta_Temp_2D, delta_N_2D, delta_2D] = RK4_2D(time_axis, Ein, T_0, N_0, a_0, solver_steps, tau_d, ...
    tau_th, tau_fc, tau_r, tau_c, eta_lin, omega_resonance_0, omega_p, m, cp, nsi, ...
    dnsi_dT, dnsi_dN, GAMMA_TH, GAMMA_FCA, V_FCA, B_TPA, c, h_bar, sigma_fca, t_sq, k_sq, zeta, phi_d, mode);


%% Import CSV

N2 = readmatrix("N_2D.csv"); 
N4 = readmatrix("N_4D.csv");
T2 = readmatrix("T_2D.csv"); 
T4 = readmatrix("T_4D.csv");
N_5 = readmatrix("N_fig_5_thomas_paper.csv"); 
T_5 = readmatrix("T_fig_5_thomas_paper.csv"); 
n2 = readmatrix("n_norm_2D.csv"); 
t2 = readmatrix("theta_norm_2D.csv");
ph_plt_fig_3 = readmatrix("phase_plot_fig_3.csv"); 

%% PHASE PLOT. Figure 3 in ""Cascadable excitability in microrings"
scatter(ph_plt_fig_3(:,1), ph_plt_fig_3(:,2)*1e23); hold on 
plot(delta_Temp, delta_N); hold on 
legend('paper', 'me')


legend({'$(\Delta T,\Delta N)(t) \; paper $', '$(\Delta T,\Delta N)(t) \; recreated$'}, ...
       'Interpreter','latex')

xlabel('$\Delta T$ (K)','Interpreter','latex')
ylabel('$\Delta N \; (m^{-3})$','Interpreter','latex')

title('$P_{in}=0.6\,\mathrm{mW},\ \delta\lambda=62\,\mathrm{pm}$', ...
      'Interpreter','latex')

grid; 



%% Time evolution. Figure 3 in "Cascadable excitability in microrings"

scatter(N2(:,1)*1e5,N2(:,2)/max(N2(:,2)) ,"filled"); hold on; 
plot(time_axis*1e6, delta_N / max(delta_N)); 

scatter(T2(:,1),T2(:,2) / max(T2(:,2)) ,"filled"); hold on; 
plot(time_axis*1e6, delta_Temp / max(delta_Temp));
xlim([0 0.5])
legend('N paper', 'N my func', 'T paper', 'T my func')
xlabel('micro second')
ylabel('Normalized values')
title(['Pin = 0.6mW,      delta lambda = 62pm.']) 
grid; 


%% Time evolution. Figure 1 in ""Simplified description of self-pulsation and excitability by thermal and free-carrier effects in semiconductor microcavities ""


plot(tau, n_norm_2D); hold on 
plot(tau, theta_norm_2D); 
scatter(n2(:,1), n2(:,2)); 
scatter(t2(:,1), t2(:,2)); 

xlim([0 10])
legend('N paper', 'N my func', 'T paper', 'T my func')
xlabel('tau')
ylabel('Normalized values')
title(['Figure 1 in paper "Simplified description of self-pulsation.... "']) 
grid; 


%% Normalized 4D equation solver 

function [E_drop, E_through, a_norm, theta_norm, n_norm, tau, a, delta_Temp, delta_N] = RK4_4D_norm(time_axis, Ein, T_0, N_0, a_0, step, tau_d, ...
    tau_th, tau_fc, tau_r, tau_c, eta_lin, omega_resonance_0, omega, m, cp, nsi, ...
    dnsi_dT, dnsi_dN, GAMMA_TH, GAMMA_FCA, GAMMA_TPA, V_FCA, V_TPA, B_TPA, c, h_bar, sigma_fca, zeta, phi_d, mode)

% Allocate arrays for memory.
theta_norm = complex(zeros(size(time_axis)));
n_norm = complex(zeros(size(time_axis)));
a_norm = complex(zeros(size(time_axis)));
a      = complex(zeros(size(time_axis)));
E_through = complex(zeros(size(time_axis))); 
E_add = complex(zeros(size(time_axis)));
E_drop = complex(zeros(size(time_axis)));


% Normalization parameters
% Make sure couling lifetime fits with the all-pass or add-drop config.
if mode == "add_drop"
    tau_ph = (1/tau_r + 2/tau_c)^(-1);
    k = tau_c/tau_r/2;
    disp(["Add-drop config assumed. tau_ph = ", num2str(tau_ph), ' and k = ' ,num2str(k), '.']); 
else
    tau_ph = (1/tau_r + 1/tau_c)^(-1);
    k = tau_c/tau_r; 
    disp(["All-pass config assumed. tau_ph = ", num2str(tau_ph), ' and k = ', num2str(k), '.']); 
end    


% Substitute 
tau_min = min(time_axis)/tau_th; 
tau_max = max(time_axis)/tau_th; 
tau_step = step/tau_th; 
tau = tau_min:tau_step:tau_max; 
tau_d_int = round(tau_d/tau_th/tau_step); % Number of time steps for delay-loop.


Q = omega_resonance_0*tau_ph/2; 
Qi = Q*(1+k)/k; 


epsilon = tau_fc/tau_th; 
mu = tau_ph/tau_th; 
delta = (omega-omega_resonance_0)*tau_ph; 
e = (1+k)/k/eta_lin; 
f = GAMMA_FCA*sigma_fca*c / (2*omega_resonance_0*abs(dnsi_dN)); 

P0_th = m*cp*nsi*((1+k)/k)^2 / (4*dnsi_dT*tau_th*eta_lin*GAMMA_TH*Qi); 
P0_el = sqrt(h_bar*omega_resonance_0^3 *nsi^3 / (4*abs(dnsi_dN)*tau_fc*GAMMA_FCA*B_TPA*c^2) ) *V_FCA/(Qi^(3/2)) * ((1+k)/k)^(3/2); 
q = P0_th / P0_el; 
 

zeta = 2*Q*dnsi_dT*tau_th*GAMMA_TH*GAMMA_TPA*B_TPA*c^2*(P0_th*tau_ph)^2 / (nsi^3 *m *cp *V_TPA);
gamma_tpa_a = GAMMA_TPA*B_TPA*c^2*P0_th*tau_ph^2 / (2000*nsi^2*V_TPA);


%Initial conditions
theta_norm(1) = 2*Q*dnsi_dT*T_0/nsi;
n_norm(1) = 2*Q*abs(dnsi_dN)*N_0/nsi; 
a_norm(1) = a_0 / sqrt(P0_th*tau_ph);


% 3D equations. d(theta)/dtau, dn/dtau and da/dtau
dtheta_dt = @(tt, theta, n, a) (-theta + abs(a).^2 * (1+e*f*n+zeta*abs(a).^2)); 
dn_dt = @(tt, theta, n, a) 1/epsilon*(-n + abs(a).^4 *q.^2); 
da_dt = @(tt, theta, n, a, p) 1/mu*(1i*(delta + theta - n) - (1+f*n + gamma_tpa_a*abs(a).^2))*a + 1i*tau_th*sqrt(p); 
h = tau_step;  % step size in tau-domain


for i = 1:length(tau)-1

    if i > tau_d_int && mode == "add_drop"
        E_add(i) = zeta*exp(-1i*phi_d)*E_through(i-tau_d_int); 
        pin = 2*abs(Ein(i)+E_add(i)).^2 / ((1+k)*P0_th*tau_ph^2); 
    else
        pin = 2*abs(Ein(i)).^2 / ((1+k)*P0_th*tau_ph^2);
    end   

    h1 = da_dt(tau(i), theta_norm(i), n_norm(i), a_norm(i), pin); 
    k1 = dn_dt(tau(i), theta_norm(i), n_norm(i), a_norm(i));
    L1 = dtheta_dt(tau(i), theta_norm(i), n_norm(i), a_norm(i));

    h2 = da_dt(tau(i) + 0.5*h, theta_norm(i)+0.5*h*L1, n_norm(i)+0.5*h*k1, a_norm(i)+0.5*h*h1, pin);
    k2 = dn_dt(tau(i) + 0.5*h, theta_norm(i) + 0.5*h*L1, n_norm(i) + 0.5*h*k1, a_norm(i)+0.5*h*h1);
    L2 = dtheta_dt(tau(i) + 0.5*h, theta_norm(i) + 0.5*h*L1, n_norm(i) + 0.5*h*k1, a_norm(i)+0.5*h*h1);

    h3 = da_dt(tau(i) + 0.5*h, theta_norm(i)+0.5*h*L2, n_norm(i)+0.5*h*k2, a_norm(i)+0.5*h*h2, pin);
    k3 = dn_dt(tau(i) + 0.5*h, theta_norm(i) + 0.5*h*L2, n_norm(i) + 0.5*h*k2, a_norm(i)+0.5*h*h2);
    L3 = dtheta_dt(tau(i) + 0.5*h, theta_norm(i) + 0.5*h*L2, n_norm(i) + 0.5*h*k2, a_norm(i)+0.5*h*h2);

    h4 = da_dt(tau(i) + h, theta_norm(i)+ h*L3, n_norm(i)+h*k3, a_norm(i)+h*h3, pin); 
    k4 = dn_dt(tau(i) + h, theta_norm(i) + h*L3, n_norm(i) + h*k3, a_norm(i)+h*h3);
    L4 = dtheta_dt(tau(i) + h, theta_norm(i) + h*L3, n_norm(i) + h*k3, a_norm(i)+h*h3);

    % Update
    a_norm(i+1)     = a_norm(i)     + (h/6)*(h1 + 2*h2 + 2*h3 + h4); 
    n_norm(i+1)     = n_norm(i)     + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
    theta_norm(i+1) = theta_norm(i) + (h/6)*(L1 + 2*L2 + 2*L3 + L4);

    a(i) = a_norm(i)*sqrt(P0_th*tau_ph);
    E_through(i) = Ein(i) + 1i*sqrt(2/tau_c)*a(i);
    E_drop(i) =  E_add(i) + 1i*sqrt(2/tau_c)*a(i);

end

% Convert to non-normalized values
delta_Temp = theta_norm* nsi/2/Q/dnsi_dT; 
delta_N = n_norm*nsi/2/Q/abs(dnsi_dN); 

end






%% Functions
function [N_mask, new_steps_per_bit] = change_N_mask(steps_per_bit,N_mask_raw)
    diff = 0;
    delta_steps = 0; 
    largest_gcd = 0; 

    while largest_gcd <= 4 && delta_steps < steps_per_bit
        
        steps_per_bit_plus = steps_per_bit + delta_steps; 
        steps_per_bit_minus = steps_per_bit - delta_steps; 
     
    while(N_mask_raw - diff*2) > 0

        N_new_plus = N_mask_raw + diff; 
        N_new_minus = N_mask_raw - diff;  

        gcd_best_plus = max(gcd(steps_per_bit_plus, N_new_plus), gcd(steps_per_bit_plus, N_new_minus));
        gcd_best_minus = max(gcd(steps_per_bit_minus, N_new_plus), gcd(steps_per_bit_minus, N_new_minus));

        
        if gcd_best_plus > largest_gcd
            largest_gcd = gcd_best_plus;
            new_steps_per_bit = steps_per_bit_plus; 

        elseif gcd_best_minus > largest_gcd 
            largest_gcd = gcd_best_minus;
            new_steps_per_bit = steps_per_bit_minus; 
        end         
        diff = diff+1;  
    end
    diff = 0; 
    delta_steps = delta_steps + 1; 
    end

    N_mask = largest_gcd; 

end








%% 2D version 

function [E_drop, E_through, a_norm, theta_norm, n_norm, tau, a, delta_Temp, delta_N, delta] = RK4_2D(time_axis, Ein, T_0, N_0, a_0, step, tau_d, ...
    tau_th, tau_fc, tau_r, tau_c, eta_lin, omega_resonance_0, omega, m, cp, nsi, ...
    dnsi_dT, dnsi_dN, GAMMA_TH, GAMMA_FCA, V_FCA, B_TPA, c, h_bar, sigma_fca, t_sq, k_sq, zeta, phi_d, mode)

% Allocate arrays for memory.
theta_norm = complex(zeros(size(time_axis)));
n_norm = complex(zeros(size(time_axis)));
delta_N = complex(zeros(size(time_axis)));
delta_Temp = complex(zeros(size(time_axis)));
a_norm = complex(zeros(size(time_axis)));
a      = complex(zeros(size(time_axis)));
E_through = complex(zeros(size(time_axis))); 
E_add = complex(zeros(size(time_axis)));
E_drop = complex(zeros(size(time_axis)));



% Substitute 
tau_min = min(time_axis)/tau_th; 
tau_max = max(time_axis)/tau_th; 
tau_step = step/tau_th; 
tau = tau_min:tau_step:tau_max; 
tau_d_int = round(tau_d/tau_th/tau_step); % Number of time steps for delay-loop.
numel(tau)


% Normalization parameters
% Make sure couling lifetime fits with the all-pass or add-drop config.
if mode == "add_drop"
    tau_ph = (1/tau_r + 2/tau_c)^(-1);
    k = tau_c/tau_r/2;
    disp(["Add-drop config assumed. tau_ph = ", num2str(tau_ph), ' and k = ' ,num2str(k), '.']); 
else
    tau_ph = (1/tau_r + 1/tau_c)^(-1);
    k = tau_c/tau_r; 
    disp(["All-pass config assumed. tau_ph = ", num2str(tau_ph), ' and k = ', num2str(k), '.']); 
end    
 
Q = omega_resonance_0*tau_ph/2;
Qi = Q*(1+k)/k;

epsilon = tau_fc/tau_th; 
delta = (omega-omega_resonance_0)*tau_ph; 
e = (1+k)/k/eta_lin; 
f = GAMMA_FCA*sigma_fca*c / (2*omega_resonance_0*abs(dnsi_dN));

P0_th = m *cp *nsi * ((1+k)/k) / (4*dnsi_dT*tau_th*eta_lin*GAMMA_TH*Q); 
P0_el = sqrt(h_bar * omega_resonance_0^3 *nsi^3  / (4*abs(dnsi_dN)*tau_fc*GAMMA_FCA*B_TPA*c^2) ) *V_FCA/(Qi^(3/2)) * ((1+k)/k)^(3/2); 
q = P0_th / P0_el;


%Initial conditions
theta_norm(1) = 2*Q*dnsi_dT*T_0/nsi;
n_norm(1) = 2*Q*abs(dnsi_dN)*N_0/nsi; 
a_norm(1) = a_0 / sqrt(P0_th*tau_ph);


% Denominator helper (same in both eqs)
% 2D equations. dn/dtau and d(theta)/dtau
denom = @(theta,n) (1 + f*n).^2 + (delta + theta - n).^2;
dn_dt = @(tt,theta,n,p) (1/epsilon) * ( -n + (p*q ./ denom(theta,n)).^2 );
dtheta_dt = @(tt,theta,n,p) ( -theta + 1*p*(1 + e*f*n) ./ denom(theta,n) );
h = tau_step;  % step size in tau-domain

for i = 1:length(tau)-1

    if i > tau_d_int && mode == "add_drop"
        E_add(i) = zeta*exp(-1i*phi_d)*E_through(i-tau_d_int); 
        pin = 2*abs(Ein(i)+E_add(i)).^2 / ((1+k)*P0_th); 
    else
        pin = 2*abs(Ein(i)).^2 / ((1+k)*P0_th);  
    end    


    k1 = dn_dt(tau(i), theta_norm(i), n_norm(i), pin);
    L1 = dtheta_dt(tau(i), theta_norm(i), n_norm(i), pin);
    k2 = dn_dt(tau(i) + 0.5*h, theta_norm(i) + 0.5*h*L1, n_norm(i) + 0.5*h*k1, pin);
    L2 = dtheta_dt(tau(i) + 0.5*h, theta_norm(i) + 0.5*h*L1, n_norm(i) + 0.5*h*k1, pin);
    k3 = dn_dt(tau(i) + 0.5*h, theta_norm(i) + 0.5*h*L2, n_norm(i) + 0.5*h*k2, pin);
    L3 = dtheta_dt(tau(i) + 0.5*h, theta_norm(i) + 0.5*h*L2, n_norm(i) + 0.5*h*k2, pin);
    k4 = dn_dt(tau(i) + h, theta_norm(i) + h*L3, n_norm(i) + h*k3, pin);
    L4 = dtheta_dt(tau(i) + h, theta_norm(i) + h*L3, n_norm(i) + h*k3, pin);

    % Update
    n_norm(i+1)     = n_norm(i)     + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
    theta_norm(i+1) = theta_norm(i) + (h/6)*(L1 + 2*L2 + 2*L3 + L4);

    a_norm(i) = 1i*sqrt(2/((1+k)*P0_th))*(Ein(i)+E_add(i)) / (1i*(delta + theta_norm(i) - n_norm(i)) + (1+f*n_norm(i)));
    a(i) = a_norm(i)*sqrt(P0_th*tau_ph);
  
    E_through(i) = Ein(i) + 1i*sqrt(2/tau_c)*a(i);
    E_drop(i) =  E_add(i) + 1i*sqrt(2/tau_c)*a(i);
end
% Convert to non-normalized values
delta_Temp = theta_norm* nsi/2/Q/dnsi_dT; 
delta_N = n_norm*nsi/2/Q/abs(dnsi_dN); 



end




