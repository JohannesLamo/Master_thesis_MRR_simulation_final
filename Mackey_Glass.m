%% Mackey-Glass generator
function [input, target] = Mackey_Glass(delay_order, total_L)

% Parameters from paper "Microring resonators with external optical
% feedback for time delay reservoir computing" by GIOVANNI DONATI.
alpha = 0.2;
beta  = 10;
gamma = 0.1;
tau_delay = 17;

dt = 0.1;              % integration step
oversampling = 3;      % as in paper


N_total = total_L + delay_order; 

% Need enough internal MG steps
N_internal = (N_total + 10) * oversampling;

delay_steps = round(tau_delay / dt);

% Initial history
x = zeros(1, N_internal + delay_steps + 1);
x(1:delay_steps+1) = 1.2;   % common MG initial condition

% Euler integration
for n = delay_steps+1 : length(x)-1

    x_tau = x(n - delay_steps);

    dx = alpha * x_tau / (1 + x_tau^beta) - gamma * x(n);

    x(n+1) = x(n) + dt * dx;
end

% Remove history and transient
x = x(delay_steps+1:end);

% Oversampling of 3:
% keep every 3rd sample
x = x(1:oversampling:end);

% Make sure length is enough
x = x(1:N_total);

% Normalize
x = (x - mean(x)) / std(x);

% Optional: rescale to [0, 0.5] if your input layer expects NARMA-like range
x_minmax = (x - min(x)) / (max(x) - min(x));
x_scaled = 0.5 * x_minmax;

% One-step-ahead prediction
input  = x_scaled(1:end-delay_order);
target = x_scaled(1+delay_order:end);

end