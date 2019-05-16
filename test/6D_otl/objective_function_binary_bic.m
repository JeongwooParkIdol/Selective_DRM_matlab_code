function binary_output = objective_function_binary_bic(x_binary, X, Y, ndv)
    variable_set = [];
    for ii = 1:ndv-2
        for jj = 1:ndv-1
            if (ii<jj) 
                variable_set = [variable_set;ii,jj];
            end
        end
    end   

% cross-term이 선택된 경우 axis term무조건 선택되도록
    criteria_cross_term = 0;
    for i3 = 1:nchoosek(ndv-1, 2)
        if( (x_binary(ndv-1+i3) == 1) && (x_binary(variable_set(i3,1))*x_binary(variable_set(i3,2)) == 0))
            criteria_cross_term = 1; 
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(criteria_cross_term == 0)
        V_data_init = [ones(size(Y,1),1)]; %consider constant value
        for i = 1:ndv-1 
            V_data_init = [V_data_init X(:,i)];
        end
     
        for i2 = 1:ndv-1
            if(x_binary(i2) == 1)
                V_data_init = [V_data_init X(:,i2).^2];
            end
        end

        for i3 = 1:nchoosek(ndv-1, 2)
            if(x_binary(ndv-1+i3) == 1)
                V_data_init = [V_data_init X(:,variable_set(i3,1)).*X(:,variable_set(i3,2))];
            end
        end

        coefficient_data_init = [regress(Y, V_data_init)];
        options = optimset('PlotFcns',@optimplotfval,'MaxFunEvals',100000, 'MaxIter', 800*(length(coefficient_data_init)+1)); %100000 ,'TolFun', 1e-3
        data = [V_data_init Y]; fun = @(x)MaxLkh_ga2(x,data);
        [visParams, MaxLkh_value] = fminsearch(fun, [coefficient_data_init;0.1*ones(1,1)'], options);
        AIC_case = 2*(MaxLkh_value + size(V_data_init,2) +1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% mod %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. AIC 
        binary_output = AIC_case
% 2. AICc 
        binary_output = AIC_case + 2*(size(V_data_init,2)+1)*(size(V_data_init,2)+2)/(size(Y,1)-size(V_data_init,2)-2)
% 3. BIC        
        binary_output = 2*(MaxLkh_value ) + (size(V_data_init,2)+1)*log(size(Y,1))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else
       binary_output = 0;
    end
end