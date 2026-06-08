function [input, target] = santa_fe(delay)
    data = load("santa_fe.mat"); 
    santa_fe = data.santaFe; 
    santa_fe_norm = (santa_fe-mean(santa_fe))/std(santa_fe); 
    input = santa_fe_norm(1:end-delay);
    target = santa_fe_norm(1+delay:end); 
end