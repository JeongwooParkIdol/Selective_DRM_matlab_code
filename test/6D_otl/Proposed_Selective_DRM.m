clc; clear all
format long e 
tic
% load('selection_dir_bic.mat')
load('selection_bic.mat')
load('DRM_results.mat')
% x = ones(1,15)
selected_uni = x(1:ndv-1); selected_bi = x(ndv:length(x));

%% DRM bivariate_mod_expectvalue
i_selected_bi = 1;
Vbi = []; Vuni1 = []; Vuni2 = [];
for i1 = 1:ndv-2
    for i2 = i1+1:ndv-1
        if(selected_bi(i_selected_bi) == 1)
            Vbi = [Vbi; [i1 i2]]; 
        end
        i_selected_bi = i_selected_bi + 1;
    end
end
for i3 = 1:ndv-1
   if(sum(sum(Vbi==i3))==0)
       Vuni1 = [Vuni1; i3];
   elseif(sum(sum(Vbi==i3))>0)
       for i4 = 1:sum(sum(Vbi==i3))-1
           Vuni2 = [Vuni2; i3];
       end
   end
end
SN = size(Vbi,1); Nd = size(Vuni2,1); Nu = size(Vuni1,1);

s_bi = struct('Vbi1',{},'Vbi2',{},'Vbi3',{},'Vbi4',{},'Vbi5',{},'Vbi6',{},'Vbi7',{},'Vbi8',{},'Vbi9',{});
g_beta_vector = V_sym;g_beta_vector(ndv) = formresults.beta1;

g_Vbi2 = s(1).(sprintf('V%d',ndv)); %performance function in V-spave(ini)
input_bi = zeros(ndv,SN); input_bi(ndv,:)= formresults.beta1;
for i = 1:SN
    T = zeros(ndv);T(Vbi(i,1),Vbi(i,1)) = 1;T(Vbi(i,2),Vbi(i,2)) = 1;T(ndv,ndv) = 1;
    g_Vbi2 = g_Vbi2 + subs(gv,V_sym, T * g_beta_vector);
    input_bi(Vbi(i,1),i) = 1;input_bi(Vbi(i,2),i) = 1;
end
for i = 1:length(Vuni1)
    g_Vbi2 = g_Vbi2 + s(1).(sprintf('V%d',Vuni1(i)));
end
for i = 1:length(Vuni2)
    g_Vbi2 = g_Vbi2 - s(1).(sprintf('V%d',Vuni2(i)));
end


% 3-point gaussian integration _ bivariate
gausian_point = 3;w_3 = [-sqrt(3) 0 sqrt(3)];weight_3 = [1/6 4/6 1/6];

input_vector =  zeros(ndv,1);
input_vector(ndv,1) = formresults.beta1;
b_1=double(subs(k,V_sym,input_vector));
beta=formresults.beta1;
pf_DRM_bi = 0;

if(SN>0)
for sub_i1 = 1:gausian_point
    input_vector_sub_d1 = zeros(ndv,1);input_vector_sub_d1(ndv) = beta;input_vector_sub_d1(1) = w_3(sub_i1);
    pf_DRM_bi_sub_d1 = normcdf(-beta+double(subs(g_Vbi2,V_sym,input_vector_sub_d1))/b_1)^(sum(Vuni2==1));
    for sub_i2 = 1:gausian_point
        input_vector_sub_d2 = zeros(ndv,1);input_vector_sub_d2(ndv) = beta;input_vector_sub_d2(2) = w_3(sub_i2);
        pf_DRM_bi_sub_d2 = normcdf(-beta+double(subs(g_Vbi2,V_sym,input_vector_sub_d2))/b_1)^(sum(Vuni2==2));
        for sub_i3 = 1:gausian_point
            input_vector_sub_d3 = zeros(ndv,1);input_vector_sub_d3(ndv) = beta;input_vector_sub_d3(3) = w_3(sub_i3);
            pf_DRM_bi_sub_d3 = normcdf(-beta+double(subs(g_Vbi2,V_sym,input_vector_sub_d3))/b_1)^(sum(Vuni2==3));
            for sub_i4 = 1:gausian_point
                input_vector_sub_d4 = zeros(ndv,1);input_vector_sub_d4(ndv) = beta;input_vector_sub_d4(4) = w_3(sub_i4);
                pf_DRM_bi_sub_d4 = normcdf(-beta+double(subs(g_Vbi2,V_sym,input_vector_sub_d4))/b_1)^(sum(Vuni2==4));
                for sub_i5 = 1:gausian_point
                    input_vector_sub_d5 = zeros(ndv,1);input_vector_sub_d5(ndv) = beta;input_vector_sub_d5(5) = w_3(sub_i5);
                    pf_DRM_bi_sub_d5 = normcdf(-beta+double(subs(g_Vbi2,V_sym,input_vector_sub_d5))/b_1)^(sum(Vuni2==5));
                  
                    input_vector =  [w_3(sub_i1); w_3(sub_i2); w_3(sub_i3); w_3(sub_i4); w_3(sub_i5); 1];
                    pf_DRM_bi_sub = 1;
                    for SN_i = 1:SN
                        pf_DRM_bi_sub = pf_DRM_bi_sub * normcdf(-beta+double(subs(g_Vbi2,V_sym,input_bi(:,SN_i).*input_vector))/b_1);
                    end
                    pf_DRM_bi = pf_DRM_bi + weight_3(sub_i1)*weight_3(sub_i2)*weight_3(sub_i3)*weight_3(sub_i4)*weight_3(sub_i5) * pf_DRM_bi_sub/pf_DRM_bi_sub_d1/ pf_DRM_bi_sub_d2/pf_DRM_bi_sub_d3/pf_DRM_bi_sub_d4/pf_DRM_bi_sub_d5;
                end
            end
        end
    end
end
else
    pf_DRM_bi = 1;
end
for i = 1:length(Vuni1)
    pf_DRM_original_sub = 0;
    for j = 1:gausian_point
        input_vector =  zeros(ndv,1);input_vector(ndv,1) = beta;
        if(selected_uni(Vuni1(i))==1)        
            input_vector(Vuni1(i)) = w_3(j);
        end
        pf_DRM_original_sub = pf_DRM_original_sub + weight_3(j)*normcdf(-beta+double(subs(g_Vbi2,V_sym,input_vector))/b_1);
    end
    pf_DRM_bi = pf_DRM_bi * pf_DRM_original_sub;
end

pf_DRM_bivariate_3point = pf_DRM_bi / normcdf(-beta)^(SN-Nd+Nu-1);
disp('PF_selectiveDRM:');
disp(pf_DRM_bivariate_3point);

toc
