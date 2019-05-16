clc; clear all
format long e 
tic

%% 1. FORM, SORM
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% mod %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%line Aver_X, Sig_X, probdata.marg-, 21-27,analysisopt.grad_flag
b1 = [50 25 0.5 1.2 0.25 50];
b2 = [150 70 3 2.5 1.2 300];
% Aver_X = [0.9337 220000 21000 0.29 24 8];
% Sig_X = [0.0459 5000 1000 0.005774 0.5 0.3] ;
Aver_X = (b1+b2)/2;
Sig_X =  (b2-b1)/6;

%Distri(1:2) = 1;
ndv = length(Aver_X); %변수의 갯수

probdata.marg(1,:) =  [ 1  Aver_X(1)   Sig_X(1)   Aver_X(1)  0 0 0 0 0];
probdata.marg(2,:) =  [ 1  Aver_X(2)   Sig_X(2)   Aver_X(2)  0 0 0 0 0];
probdata.marg(3,:) =  [ 1  Aver_X(3)   Sig_X(3)   Aver_X(3)  0 0 0 0 0];
probdata.marg(4,:) =  [ 1  Aver_X(4)   Sig_X(4)   Aver_X(4)  0 0 0 0 0];
probdata.marg(5,:) =  [ 1  Aver_X(5)   Sig_X(5)   Aver_X(5)  0 0 0 0 0];
probdata.marg(6,:) =  [ 1  Aver_X(6)   Sig_X(6)   Aver_X(6)  0 0 0 0 0];

probdata.correlation = eye(ndv); %4*4행렬I만듬

probdata.parameter = distribution_parameter(probdata.marg);

%eq_170816
% gfundata(1).expression = ' ((3*x(1)*x(2)*385.82*(x(5)-x(6))/(x(4)*(x(3)*2*pi/60)^2*(x(5)^3-x(6)^3)))^0.5 - 0.37473)'; % -( 여기다가 함수 써야됨)

gfundata(1).expression = ' ((12*x(2)/(x(1)+x(2))+0.74)*x(6)*(x(5)+9)/(x(6)*(x(5)+9)+x(3)) + 11.35*x(3)/(x(6)*(x(5)+9)+x(3)) + 0.74*x(3)*x(6)*(x(5)+9)/x(4)/(x(6)*(x(5)+9)+x(3))-4)'; % -( 여기다가 함수 써야됨)
gfundata(1).evaluator = 'basic';
gfundata(1).type = 'expression'; 
gfundata(1).parameter = 'no';
gfundata(1).dgdq={'-(3*x(1)^2 + x(2) + x(4))';
                  '-(2*x(2) + x(1) - 2*x(2)*x(3)*x(4))';
                  '-(-1-x(4)-x(2)^2*x(4))';
                  '(1 + x(1) - x(3) - x(2)^2 * x(3))'}    ;   

% FORM analysis options8613083
analysisopt.ig_max    = 100;
analysisopt.il_max    = 5;
analysisopt.i_max = 100; % Maximum number of iterations allowed in the search algorithm
analysisopt.e1 = 0.001; % Tolerance on how close design point is to limit-state surface
analysisopt.e2 = 0.001; % Tolerance on how accurately the gradient points towards
                        % the origin
analysisopt.step_code = 0; % 0: step size by Armijo rule
                            % otherwise: given value is the step size
analysisopt.Recorded_u = 1; % 0: u-vector not recorded at all iterations
                            % 1: u-vector recorded at all iterations
analysisopt.Recorded_x = 1; % 0: x-vector not recorded at all iterations
                            % 1: x-vector recorded at all iterations
analysisopt.grad_flag = 'ffd';

% IC analysis options
analysisopt.rand_generator = 1; % 0: default rand matlab function
                                % 1: Mersenne Twister (to be preferred)
analysisopt.num_sim = 1e8; % Number of samples

analysisopt.stdv_sim  = 1;

analysisopt.target_cov = 0.001;
analysisopt.lowRAM = 0; % 1: memory savings allowed
                        % 0: no memory savings allowed

femodel = 0;
randomfield.mesh = 0;

% performance function

% IC : dspt, CMCS : origin
analysisopt.sim_point = 'dspt' ;
% the nubmer of sampling points
%analysisopt.num_sim   = 1e04;

%% FORM
[formresults] = form(1,probdata,analysisopt,gfundata,femodel,randomfield)
Grad_Total = formresults.grad_G;
X_Total = formresults.u;

%% Conventional SORM
sormresults = sorm_curvature_fitting(1,formresults,probdata,analysisopt,gfundata,femodel,randomfield)

% formresults.alpha
%%
%simulationresults = simulation(1,probdata,analysisopt,gfundata,femodel,randomfield)
save('FORM_results.mat')
disp('PF_FORM:');
disp(formresults.pf1);