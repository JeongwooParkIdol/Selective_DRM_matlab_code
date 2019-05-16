%% AIC, AICc, BIC
clc; clear all
format long e 
tic

load('FORM_results.mat')
load('DRM_results.mat')

num_sample = size((formresults.u_eval),2);

X = (R\formresults.u_eval)';
Y = [];
for y_k=1:num_sample
    Y=[Y;eval(subs(gv,V_sym,X(y_k,:)')-b1*(X(y_k,ndv)-formresults.beta1))];
end

num_aixs = ndv-1;
num_cross = nchoosek(ndv-1, 2);
num_binary = num_aixs + num_cross;
IntCon=[1:num_binary];

%Objective function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% mod %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vfun=@(x)objective_function_binary_bic(x, X, Y, ndv);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Initial population
ini=zeros(1,num_binary);

%lb,ub
lb=ones(1,num_binary)*(-0.5);
ub=ones(1,num_binary)*(+1.5);

%GA option

options = gaoptimset('InitialPopulation',ini,'PopulationSize',100,'Generations',200,'PlotFcns',@gaplotbestfun,'UseParallel','always');

%GA
[x,fval,exitflag,output] = ga(vfun,num_binary,[],[],[],[],lb,ub,'nonlcon',IntCon, options)

save('selection_dir_bic.mat')