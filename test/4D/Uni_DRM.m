clc; clear all
format long e 
tic

%% MPP-based DRM
load('FORM_results.mat')
for i = 1:15
    syms (sprintf('U%d',i));
    syms (sprintf('u%d',i));
    syms (sprintf('V%d',i));
    syms (sprintf('x%d',i));
end
U_vector_sym = [U1; U2; U3; U4; U5; U6; U7; U8; U9; U10; U11; U12; U13; U14; U15]; U_sym = U_vector_sym(1:ndv);
V_vector_sym = [V1; V2; V3; V4; V5; V6; V7; V8; V9; V10; V11; V12; V13; V14; V15]; V_sym = V_vector_sym(1:ndv); 
x_vector_sym = [x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15]; x_sym = x_vector_sym(1:ndv); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% mod %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g = ((4*(5*U4+100).^3./(2.8*3.3*(1450000*U1+29000000)).*sqrt(((100*U3+1000)/3.3^2).^2+((100*U2+500)/2.8^2).^2)-2.5));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R = gram_schmidt(formresults.alpha)';
U=R*V_sym;
gv=subs(g,U_sym,U);
k=diff(gv,V_sym(ndv));
input_vector =  zeros(ndv,1);
input_vector(ndv,1) = formresults.beta1;
b1=double(subs(k,V_sym,input_vector));


%% DRM univariate
s = struct('V1',{},'V2',{},'V3',{},'V4',{},'V5',{},'V6',{},'V7',{},'V8',{},'V9',{});
g_beta_vector = V_sym;
mpp = formresults.alpha * formresults.beta1; %% bivariate효과 알아보기 위해 mpp로 계산
g_beta_vector(ndv) = formresults.beta1;
g_V = 0; %V-space 에서의 performance function
for i = 1:ndv
    
    T = zeros(ndv);
    T(i,i) = 1;
    T(ndv,ndv) = 1;
    g_beta_vector_mod = T * g_beta_vector;
    if(i == ndv)
        s(1).(sprintf('V%d',i))=subs(k,V_sym,g_beta_vector_mod)*(V_sym(ndv)-formresults.beta1);
    else
        
        s(1).(sprintf('V%d',i)) = subs(gv,V_sym, g_beta_vector_mod);
    end
    g_V = g_V + s(1).(sprintf('V%d',i));
end
g_V_mppbased = g_V;

%3-point
gausian_point = 3;
sum_matrix = zeros(ndv-1,gausian_point);
w_3 = [-sqrt(3) 0 sqrt(3)];
weight_3 = [0.166667 0.666667 0.166667];

input_vector =  zeros(ndv,1);
input_vector(ndv,1) = formresults.beta1;
b1=double(subs(k,V_sym,input_vector));
b=formresults.beta1;
pf_DRM_original = 1;

for i = 1:ndv-1
     pf_DRM_original_sub = 0;
    for j = 1:gausian_point
        input_vector =  zeros(ndv,1);
        input_vector(i,1) = w_3(j);
        input_vector(ndv,1) = formresults.beta1;
        sum_matrix(i,j) = double(subs(g_V,V_sym,input_vector))  ;
        pf_DRM_original_sub = pf_DRM_original_sub + weight_3(j)*normcdf(-b+sum_matrix(i,j)/b1);
    end

    pf_DRM_original = pf_DRM_original * pf_DRM_original_sub;
end

pf_DRM_original_3point = pf_DRM_original / normcdf(-b)^(ndv-2);

%7-point
gausian_point = 7;
sum_matrix = zeros(ndv-1,gausian_point);
w_7=[-3.750439717725743  -2.366759410734542  -1.154405394739968     0    1.154405394739968   2.366759410734542   3.750439717725743];
weight_7=[0.000548268855972   0.030757123967586   0.240123178605013   0.457142857142857   0.240123178605013   0.030757123967586   0.000548268855972];

input_vector =  zeros(ndv,1);
input_vector(ndv,1) = formresults.beta1;
b1=double(subs(k,V_sym,input_vector));
b=formresults.beta1;
pf_DRM_original_7point = 1;

for i = 1:ndv-1
    pf_DRM_original_sub = 0;
    for j = 1:gausian_point
        input_vector =  zeros(ndv,1);
        input_vector(i,1) = w_7(j);
        input_vector(ndv,1) = formresults.beta1;
        sum_matrix(i,j) = double(subs(g_V,V_sym,input_vector));
        pf_DRM_original_sub = pf_DRM_original_sub + weight_7(j)*normcdf(-b+sum_matrix(i,j)/b1);
    end
    pf_DRM_original_7point = pf_DRM_original_7point * pf_DRM_original_sub;
end
pf_DRM_original_7point = pf_DRM_original_7point / normcdf(-b)^(ndv-2);

disp('PF_uniDRM_3pts:');
disp(pf_DRM_original_3point);
disp('PF_uniDRM_7pts:');
disp(pf_DRM_original_7point);
save('DRM_results.mat')