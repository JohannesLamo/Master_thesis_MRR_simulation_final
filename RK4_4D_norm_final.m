function [E_drop, a, delta_lambda_T_plus_FC] = RK4_4D_norm_final(time_axis, Ein, T_0, N_0, a_0, step, tau_d, ...
    tau_th, tau_fc, tau_r, tau_c, eta_lin, omega_resonance_0, omega, m, cp, nsi, ...
    dnsi_dT, dnsi_dN, GAMMA_TH, GAMMA_FCA, GAMMA_TPA, V_FCA, V_TPA, B_TPA, c, h_bar, sigma_fca, zeta, phi_d, N_mask, mode)

% This function solves the 3D system of differential equations for an MRR.
% Below are a list of the inputs and outputs of this function 

% ========================= FUNCTION INPUTS =========================

% time_axis        : [vector] Physical time axis [s]
% Ein              : [complex vector] Input optical field (pump + modulation) [sqrt(W)]
% T_0              : [scalar] dT at time = 0 [K]
% N_0              : [scalar] dN at time = 0 [m^(-3)]
% a_0              : [scalar] complex modal amplitude inside MRR at time = 0 [sqrt(J)]
% step             : [scalar] Solver time step in physical domain [s]
% tau_d            : [scalar] Delay time in delay loop [s]
%
% tau_th           : [scalar] Thermal relaxation time [s]
% tau_fc           : [scalar] Free carrier lifetime [s]
% tau_r            : [scalar] Intrinsic loss lifetime [s]
% tau_c            : [scalar] Coupling lifetime [s]
%
% eta_lin          : [scalar] Fraction of linear loss due to absorption [~]
% omega_resonance_0: [scalar] Cold cavity resonance angular frequency [rad/s]
% omega            : [scalar] Pump angular frequency [rad/s]
%
% m                : [scalar] Effective mass of MRR [kg]
% cp               : [scalar] Specific heat capacity [J/(kg·K)]
% nsi              : [scalar] Refractive index of silicon [~]
%
% dnsi_dT          : [scalar] Thermo-optic coefficient [K^-1]
% dnsi_dN          : [scalar] Free-carrier dispersion coefficient [m^3]
%
% GAMMA_TH         : [scalar] Thermal confinement factor [~]
% GAMMA_FCA        : [scalar] Free-carrier absorption confinement factor [~]
% GAMMA_TPA        : [scalar] TPA confinement factor [~]
%
% V_FCA            : [scalar] Effective mode volume for FCA [m^3]
% V_TPA            : [scalar] Effective mode volume for TPA [m^3]
% B_TPA            : [scalar] TPA coefficient [m/W]
%
% c                : [scalar] Speed of light [m/s]
% h_bar            : [scalar] Reduced Planck constant [J·s]
% sigma_fca        : [scalar] FCA cross-section [m^2]
%
% zeta             : [scalar] coupling factor into delay loop [~]
% phi_d            : [scalar] Phase shift in delay loop [rad]
%
% mode             : [string] "add_drop" or "all_pass" configuration


% ========================= FUNCTION OUTPUTS =========================
%
% E_drop                    : [complex vector] Optical field at drop port [sqrt(W)]

% a                         : [complex vector] Physical complex modal amplitude [sqrt(J)]
% delta_lambda_T_plus_FC    : [vector] Change in resonance wavelength [m]
% ===================================================================


% Allocate complex arrays
theta_norm = complex(zeros(size(time_axis)));
n_norm = complex(zeros(size(time_axis)));
a_norm = complex(zeros(size(time_axis)));
a      = complex(zeros(size(time_axis)));
E_through = complex(zeros(size(time_axis))); 
E_add = complex(zeros(size(time_axis)));
E_drop = complex(zeros(size(time_axis))); 
delta_lambda_T_plus_FC = zeros(size(time_axis)); 

% Convert the physical time axis into the normalized time variable tau.
% tau_step is the RK4 step size in the normalized time domain.
tau_min = min(time_axis)/tau_th; 
tau_max = max(time_axis)/tau_th; 
tau_step = step/tau_th; % step size in tau-domain
tau = tau_min:tau_step:tau_max; 
h = tau_step;  % step size in tau-domain

% Convert the physical delay time into an integer number of RK4 steps.
tau_d_int = round(tau_d/tau_th/tau_step); % Number of time steps for delay-loop.

% Photon lifetime in the cavity. Uneffected by mode being = "add_drop" or "all_pass"
tau_ph = (1/tau_r + 1/tau_c)^(-1);

% Ratio between coupling lifetime and intrinsic loss lifetime. Uneffected by mode being = "add_drop" or "all_pass"
k = tau_c/tau_r; 

Q = omega_resonance_0*tau_ph/2;                                             %Loaded quality factor
Qi = Q*(1+k)/k;                                                             %Intrinsic quality factor

epsilon = tau_fc/tau_th;                                                    %Ratio between free-carrier and thermal time scales
mu = tau_ph/tau_th;                                                         %Ratio between photon and thermal time scales
delta = (omega-omega_resonance_0)*tau_ph;                                   %Normalized detuning frequnecy
e = (1+k)/k/eta_lin;                                                        %Related to FCA-induced heating
f = GAMMA_FCA*sigma_fca*c / (2*omega_resonance_0*abs(dnsi_dN));             %Ratio of FCA to FCD 


P0_el = sqrt(h_bar*omega*V_FCA^2*nsi^3 / (Q * abs(dnsi_dN) *GAMMA_FCA *tau_fc*c^2*B_TPA*tau_ph^2));
P0_th = m*cp*nsi / (4*Q*abs(dnsi_dT)*tau_th*GAMMA_TH*eta_lin/tau_r *tau_ph); 

% Ratio between the thermal and electrical normalization powers.
q = P0_th / P0_el; 
 
  
eta = 2*Q*dnsi_dT*tau_th*GAMMA_TH*GAMMA_TPA*B_TPA*c^2*(P0_th*tau_ph)^2 / (m*cp*nsi^3*V_TPA);     %Related to TPA-induced heating. Can be found diff.eq. for theta
gamma_tpa = GAMMA_TPA*B_TPA*c^2*P0_th*tau_ph^2 / (2*nsi^2*V_TPA);                                %Normalized TPA loss term in field equation. Can be found in diff.eq for a

% Initial conditions in normalized variables.
theta_norm(1) = 2*Q*dnsi_dT*T_0/nsi;
n_norm(1) = 2*Q*abs(dnsi_dN)*N_0/nsi; 
a_norm(1) = a_0 / sqrt(P0_th*tau_ph);

dtheta_dt = @(tt, theta, n, a) -theta + abs(a).^2*(1+e*f*n + eta*abs(a).^2);  
dn_dt = @(tt, theta, n, a) 1/epsilon*(-n  + abs(a).^4 * q^2 ); 

if mode == "add_drop"
    x = (tau_c + 2*tau_r) / (tau_c + tau_r); 
    %disp("Add-drop config assumed. "); 
    da_dt = @(tt, theta, n, a, ein, eadd) 1/mu*(1i*(delta + theta -n)   - (x + f*n + gamma_tpa*abs(a).^2)).*a ....
    + 1i* tau_th*  sqrt(2/tau_c)*1/sqrt(P0_th*tau_ph)*(ein+eadd);
else
    x = 1;  
    %disp("All-pass config assumed."); 
    da_dt = @(tt, theta, n, a, ein, eadd) 1/mu*(1i*(delta + theta -n) - (x + f*n + gamma_tpa*abs(a).^2))*a ....
    + 1i* tau_th* sqrt(2/tau_c)*1/sqrt(P0_th*tau_ph)*ein;
end

% Main RK4 loop over all normalized time steps.
for i = 1:length(tau)-1

    % If the current time index is beyond the delay length and the mode is
    % add-drop, include the delayed feedback field from the through port.
    if i > tau_d_int 
        E_add(i) = sqrt(zeta)*exp(-1i*phi_d)* E_through(i-tau_d_int); 
    else
        E_add(i) = 0; 
    end 

    % Compute the four RK4 stages for the intracavity field, carrier density,
    % and temperature equations.
    h1 = da_dt(tau(i), theta_norm(i), n_norm(i), a_norm(i), Ein(i), E_add(i)); 
    k1 = dn_dt(tau(i), theta_norm(i), n_norm(i), a_norm(i));
    L1 = dtheta_dt(tau(i), theta_norm(i), n_norm(i), a_norm(i));

    h2 = da_dt(tau(i) + 0.5*h, theta_norm(i)+0.5*h*L1, n_norm(i)+0.5*h*k1, a_norm(i)+0.5*h*h1, Ein(i), E_add(i));
    k2 = dn_dt(tau(i) + 0.5*h, theta_norm(i) + 0.5*h*L1, n_norm(i) + 0.5*h*k1, a_norm(i)+0.5*h*h1);
    L2 = dtheta_dt(tau(i) + 0.5*h, theta_norm(i) + 0.5*h*L1, n_norm(i) + 0.5*h*k1, a_norm(i)+0.5*h*h1);

    h3 = da_dt(tau(i) + 0.5*h, theta_norm(i)+0.5*h*L2, n_norm(i)+0.5*h*k2, a_norm(i)+0.5*h*h2, Ein(i), E_add(i));
    k3 = dn_dt(tau(i) + 0.5*h, theta_norm(i) + 0.5*h*L2, n_norm(i) + 0.5*h*k2, a_norm(i)+0.5*h*h2);
    L3 = dtheta_dt(tau(i) + 0.5*h, theta_norm(i) + 0.5*h*L2, n_norm(i) + 0.5*h*k2, a_norm(i)+0.5*h*h2);

    h4 = da_dt(tau(i) + h, theta_norm(i)+ h*L3, n_norm(i)+h*k3, a_norm(i)+h*h3, Ein(i), E_add(i)); 
    k4 = dn_dt(tau(i) + h, theta_norm(i) + h*L3, n_norm(i) + h*k3, a_norm(i)+h*h3);
    L4 = dtheta_dt(tau(i) + h, theta_norm(i) + h*L3, n_norm(i) + h*k3, a_norm(i)+h*h3);

    % Update the normalized intracavity field, carrier density, and temperature using RK4.
    a_norm(i+1)     = (a_norm(i)    + (h/6)*(h1 + 2*h2 + 2*h3 + h4)); 
    n_norm(i+1)     = n_norm(i)     + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
    theta_norm(i+1) = theta_norm(i) + (h/6)*(L1 + 2*L2 + 2*L3 + L4);

    % Convert the normalized intracavity field back to physical units.
    a(i) = a_norm(i)*sqrt(P0_th*tau_ph);

    % Compute the through-port output field.
    E_through(i) = Ein(i) +  1i*sqrt(2/tau_c)*a(i);

    % Compute the drop-port output field.
    E_drop(i) =  E_add(i) +  1i*sqrt(2/tau_c)*a(i); 

    omega_T_plus_FC = omega_resonance_0 - theta_norm(i)/tau_ph + n_norm(i)/tau_ph;

    delta_lambda_T_plus_FC(i) = calc_delta_lambda_final(omega_T_plus_FC, omega_resonance_0, c); 
end

a(end) = a_norm(end)*sqrt(P0_th*tau_ph);
E_drop(end) = E_add(end) + 1i*sqrt(2/tau_c)*a(end);

omega_T_plus_FC = omega_resonance_0 - theta_norm(end)/tau_ph + n_norm(end)/tau_ph;

delta_lambda_T_plus_FC(end) = calc_delta_lambda_final(omega_T_plus_FC, omega_resonance_0, c);

% Convert the normalized temperature and carrier density back to physical units.
%delta_Temp = theta_norm* nsi/2/Q/dnsi_dT; 
%delta_N = n_norm*nsi/2/Q/abs(dnsi_dN); 

end