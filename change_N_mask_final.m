function [N_mask, new_steps_per_bit] = change_N_mask_final(steps_per_bit,N_mask_raw)

    diff = 0;                 % Offset used to modify N_mask (both + and - directions)
    delta_steps = 0;          % Offset used to modify steps_per_bit
    largest_gcd = 0;          % Stores the best (largest) GCD found

    % Check if steps_per_bit is already divisible by N_mask_raw
    if mod(steps_per_bit,N_mask_raw) == 0 && steps_per_bit > 0 && N_mask_raw > 0
        N_mask = N_mask_raw; 
        new_steps_per_bit = steps_per_bit; 
        return                % Exit function early if condition is satisfied
    end

    % Loop over increasing adjustments of steps_per_bit
    while largest_gcd <= 4 && delta_steps < steps_per_bit
        
        steps_per_bit_plus = steps_per_bit + delta_steps;   % Try increasing steps_per_bit
        steps_per_bit_minus = steps_per_bit - delta_steps;  % Try decreasing steps_per_bit
     
        % Loop over possible adjustments of N_mask
        while (N_mask_raw - diff*2) > 0

            N_new_plus = N_mask_raw + diff;    % Increase N_mask
            N_new_minus = N_mask_raw - diff;   % Decrease N_mask  

            % Compute best GCD for increased steps_per_bit
            gcd_best_plus = max(gcd(steps_per_bit_plus, N_new_plus), ...
                                gcd(steps_per_bit_plus, N_new_minus));

            % Compute best GCD for decreased steps_per_bit
            gcd_best_minus = max(gcd(steps_per_bit_minus, N_new_plus), ...
                                 gcd(steps_per_bit_minus, N_new_minus));

            % Update best solution if larger GCD is found (increase case)
            if gcd_best_plus > largest_gcd
                largest_gcd = gcd_best_plus;
                new_steps_per_bit = steps_per_bit_plus; 

            % Update best solution if larger GCD is found (decrease case)
            elseif gcd_best_minus > largest_gcd 
                largest_gcd = gcd_best_minus;
                new_steps_per_bit = steps_per_bit_minus; 
            end         

            diff = diff + 1;   % Increase N_mask search range
        end

        diff = 0;                    % Reset N_mask offset
        delta_steps = delta_steps + 1; % Increase steps_per_bit offset
    end

    % Set final N_mask as the best GCD found
    N_mask = largest_gcd; 

end
