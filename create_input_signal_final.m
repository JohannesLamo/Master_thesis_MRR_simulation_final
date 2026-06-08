function [Ein, input, target] = create_input_signal_final(total_L, steps_per_bit, N_mask, Pin_avg_lin, bias, type, seed, P, delay_santa_fe, delay_mackey_glass_order)               
    
    % This function does the following: 
    % - Creates an input signal 
    % - Creates a target function
    % - converts the input signal to an electric field envelope 
    % which can be sent into the reservoir layer 

   
    % INPUTS:
    %   total_L                   - Total number of input symbols.
    %
    %   steps_per_bit             - Number of simulation time steps used to
    %                               represent a single input symbol.
    %
    %   N_mask                    - Number of virtual nodes. The mask is
    %                               divided into N_mask segments.
    %
    %   Pin_avg_lin               - Desired average optical input power [W]
    %                               after masking and biasing.
    %
    %   bias                      - DC bias added to the masked input
    %                               sequence. Used to ensure a sufficiently
    %                               large optical carrier and small relative
    %                               modulation depth.
    %
    %   type                      - String specifying the benchmark task:
    %                               "ones", "NARMA_P",
    %                               "santa_fe", or "mackey_glass".
    %
    %   seed                      - Random seed used when generating NARMA
    %                               input sequences.
    %
    %   P                         - Order of the NARMA benchmark.
    %                               Only used when type = "NARMA_P".
    %
    %   delay_santa_fe            - Prediction delay used for the Santa Fe
    %                               benchmark.
    %
    %   delay_mackey_glass_order  - Prediction order/delay used for the
    %                               Mackey-Glass benchmark.
    %
    %
    % OUTPUTS:
    %   Ein                       - Optical field envelope sent into the
    %                               reservoir. Normalized according to
    %                               |Ein|^2 = Pin.
    %
    %   input                     - Unmasked input sequence u(n).
    %
    %   target                    - Desired target sequence y(n) used for
    %                               training and testing the readout layer.
        
    if type == "ones"
        input = ones(1,total_L);                    
        target = ones(1,total_L);                                      
    elseif type == "NARMA_P"
        [input, target] = NARMA_P_final(total_L,P,seed);  
    elseif type == "santa_fe"
        [input, target] = santa_fe(delay_santa_fe); 
    elseif type == "mackey_glass"
        [input, target] = Mackey_Glass(delay_mackey_glass_order, total_L); 
    else
        error("No valid data type.")
    end    


    % After the input and target signals are defined, we mask the input signal
    mask_seed = 343;                                         % Masking seed
    rng(mask_seed);                                          % Used for reproducability 
    mask = rand(1,N_mask);                                   % Create mask of dimensions 1 X N_mask 
    M = steps_per_bit/N_mask;                                % Number of time steps per virtual node
    mask_one_bit = repelem(mask, M);                         % Replicate the mask so that we have M steps per virtual node


    masked_data_sequence = zeros(1, total_L*steps_per_bit);  % Allocate for memory

    % Multiply bit sequence by mask
    for k = 1:total_L                                        
        idx_start = (k-1)*steps_per_bit + 1;
        idx_end   = k*steps_per_bit;  
        masked_data_sequence(idx_start:idx_end) = input(k)*mask_one_bit;
    end

    data_sequence = masked_data_sequence + bias;             % Add bias to masked_data.
    
    Pin = data_sequence/mean(data_sequence)*Pin_avg_lin;     % Assume linear modulator. We want the average intensity = Pin_avg_lin.
    Ein = sqrt(Pin);                                         % Electric field envelope encoded with input data. 
end



