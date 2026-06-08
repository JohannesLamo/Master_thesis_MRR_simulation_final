function NMSE = calc_NMSE_final(y, y_hat)
    %Remove start and end in case either of these values are implemented
    %incorrectly
    y = y(2:end-1); 
    y_hat = y_hat(2:end-1);

    NMSE = mean((y-y_hat).^2)/var(y);
   
    if NMSE > 1
        NMSE = 1;
    elseif isnan(NMSE)
        NMSE = 1; 
    end    
end