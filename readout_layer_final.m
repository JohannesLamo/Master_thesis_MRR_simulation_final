function [y_train, y_test, y_train_hat, y_test_hat, input_vec, input_train, ...
    input_test, target_vec, X_train_test, time_axis_reduced_train_test, delta_lambda_T_plus_FC_reduced, time_axis_reduced_test] = readout_layer_final(E_drop, target, ...
    time_axis, steps_per_bit, N_mask, Warm1, Warm2, Ltrain, Ltest, total_L, input, LAMBDA, delta_lambda_T_plus_FC)

    % READOUT_LAYER_FINAL Processes the reservoir output and trains/tests
    % a linear readout layer using ridge regression.

    % INPUTS:
    %   E_drop                  - Optical field envelope at the drop port.
    %
    %   target                  - Desired target sequence used for training
    %                             and testing.
    %
    %   time_axis               - Full simulation time axis.
    %
    %   steps_per_bit           - Number of solver time steps per input symbol.
    %
    %   N_mask                  - Number of virtual nodes per input symbol.
    %
    %   Warm1                   - Number of initial symbols discarded before
    %                             training.
    %
    %   Warm2                   - Number of symbols discarded between training
    %                             and testing.
    %
    %   Ltrain                  - Number of symbols used for training.
    %
    %   Ltest                   - Number of symbols used for testing.
    %
    %   total_L                 - Total number of input symbols.
    %
    %   input                   - Original unmasked input sequence.
    %
    %   LAMBDA                  - Ridge regression regularization parameter.
    %
    %   delta_lambda_T_plus_FC  - Resonance wavelength shift caused by thermal
    %                             and free-carrier effects.
    %
    % OUTPUTS:
    %   y_train                 - Target values used during training.
    %
    %   y_test                  - Target values used during testing.
    %
    %   y_train_hat             - Predicted training output from the readout.
    %
    %   y_test_hat              - Predicted testing output from the readout.
    %
    %   input_vec               - Combined input sequence for training and
    %                             testing.
    %
    %   input_train             - Input sequence corresponding to the training
    %                             interval.
    %
    %   input_test              - Input sequence corresponding to the testing
    %                             interval.
    %
    %   target_vec              - Combined target sequence for training and
    %                             testing.
    %
    %   X_train_test            - Combined reservoir state matrix for training
    %                             and testing.
    %
    %   time_axis_reduced_train_test
    %                           - Reduced time axis corresponding to
    %                             X_train_test.
    %
    %   delta_lambda_T_plus_FC_reduced
    %                           - Reduced resonance wavelength shift for the
    %                             testing interval.
    %
    %   time_axis_reduced_test  - Reduced time axis corresponding to X_test.


    % Make sure that the number of time steps per bit can be evenly divided
    % into N_mask virtual nodes.
    assert(mod(steps_per_bit, N_mask) == 0, 'steps_per_bit must be divisible by N_mask');

    % Number of solver steps per virtual node.
    M = steps_per_bit / N_mask;

    % Convert the optical field at the drop port into optical intensity.
    P_drop = abs(E_drop).^2;


    % Make sure that the total number of power samples can be evenly grouped
    % into segments of length M.
    assert(mod(numel(P_drop), M) == 0, 'numel(P_drop) must be divisible by M');

    % Make sure that the total number of time samples can be evenly grouped
    % into segments of length M.
    assert(mod(numel(time_axis), M) == 0, 'numel(time_axis) must be divisible by M');

    % Reduce the original time axis by averaging over M samples.
    time_axis_reduced = mean(reshape(time_axis, M, []), 1);
    delta_lambda_T_plus_FC_reduced = mean(reshape(delta_lambda_T_plus_FC, M, []), 1);

    % Reduce the vector size of drop-port intensity by averaging over M samples.
    X_drop = mean(reshape(P_drop, M, []), 1);

    % Reshape the averaged optical intensity values into a matrix of size
    % [total number of bits] x [number of virtual nodes].
    X_drop = reshape(X_drop, N_mask, total_L).';
   
    % Extract the training feature matrix after removing the initial warm-up.
    X_train = X_drop(Warm1+1:Warm1+Ltrain,:); 

    % Extract the testing feature matrix after the training set and the
    % second warm-up interval.
    X_test = X_drop(Warm1+Ltrain+Warm2+1:end,:);
 
    % Extract the training target values corresponding to X_train.
    y_train = target(Warm1+1 : Warm1+Ltrain);

    % Make sure the training target is a column vector.
    y_train = y_train(:);

    % Extract the testing target values corresponding to X_test.
    y_test = target(Warm1+Ltrain+Warm2+1 : Warm1+Ltrain+Warm2+Ltest);

    % Make sure the testing target is a column vector.
    y_test = y_test(:);

    % Train the linear readout weights using ridge regression.
    % The first element of W is the bias term.
    W = ridge(y_train, X_train, LAMBDA, 0);

    % Compute the predicted training output using the trained readout weights.
    y_train_hat = W(1) + X_train * W(2:end);

    % Compute the predicted testing output using the trained readout weights.
    y_test_hat = W(1) + X_test * W(2:end);

    % Replace the final predicted test value with the previous one, because
    % the last predicted value is quite off. 
    y_test_hat(end) = y_test_hat(end-1); 

    
    % Save into file
    input_train = input(Warm1+1 : Warm1+Ltrain);
    input_test = input(Warm1+Ltrain+Warm2+1 : Warm1+Ltrain+Warm2+Ltest); 
    
    input_train = input_train(:);
    input_test = input_test(:);
    
    input_vec = [input_train;input_test]; 
    input_vec = input_vec(:);

    target_vec = [y_train;y_test]; 
    target_vec = target_vec(:); 

    X_train_test = [X_train; X_test]; 

    time_axis_reduced_train_test = time_axis_reduced(1:numel(X_train_test));
    time_axis_reduced_train_test = time_axis_reduced_train_test(:); 

    time_axis_reduced_test = time_axis_reduced(1:numel(X_test));
    time_axis_reduced_test = time_axis_reduced_test(:);

    % Only store reduced wavelength data for test set. 
    delta_lambda_T_plus_FC_reduced = delta_lambda_T_plus_FC_reduced(Warm1+Ltrain+Warm2+1:end);
end



