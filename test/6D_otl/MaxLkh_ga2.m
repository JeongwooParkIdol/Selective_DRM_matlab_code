function [F] = MaxLkh_ga2(params, data)

%Compute model values for data
% params: parameters(sigma  
% data    
    F = 0;


    for i = 1:size(data,1) %sample°³¼ö
%        F = F * 1/(2*pi*params(length(params))^2)^0.5 * exp(-1/2/(params(length(params)))^2 * (data(i,size(data,2)) - data(i,1:size(data,2)-1) * params(1:length(params)-1))^2);
        F = F + 1/2*log(2*pi*params(length(params))^2) + 1/(2*params(length(params))^2) * (data(i,size(data,2)) - data(i,1:size(data,2)-1) * params(1:length(params)-1))^2;
    end

end

