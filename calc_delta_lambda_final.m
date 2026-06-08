% This function calulates the current change in resonance wavelength
function delta_lambda = calc_delta_lambda_final(omega, omega_resonance_0, c)
    % INPUTS: 
    % -omega: pump frequency
    % -omega_resonance_0: Cold cavity resonance
    % -c: Speed of light in vacuum

    % OUTPUT: 
    % -delta_lambda: Change in resonance wavelength
    lambda = 2*pi*c/omega; 
    lambda_resonance = 2*pi*c/omega_resonance_0;
    delta_lambda = lambda_resonance - lambda; 
end