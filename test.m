clear all; close all;
addpath('Roomsimove_1.4/Roomsimove/')

% [h,fe] = audioread('data/C4DM_RIR_Omni/00x00y.wav');
% h = roomsimove_single('room_sensor_config.txt',[1;1.5;1]);
% h = load('h.mat');
% h = h.h;
fe = 16000;
rt60 = 0.1;
% abs = RT60toA([4.45;3.55;2.5], rt60)/fe;
lambda = 4*pi*200*340^3/(fe^3);
absorp = 3*log(10)/(rt60*fe);
var_h = 1e-8;
L_g = 200;
L_h = 2000;
L_h_true = L_h - L_g + 1;
beta = (1:L_h_true)';
alpha = (lambda*(1:L_h_true).^2)';

[h,g] = synthesiseH(L_h, L_g, absorp, lambda, var_h);

algo = VEM('h', h, 'fe', fe, 'var_h', var_h, 'L_g', L_g, 'lambda', lambda, 'abs', absorp, 'beta', beta);%, 'alpha', alpha,  'beta', beta);%, 'lambda', lambda, 'abs', abs, 'alpha', alpha);%, 'lambda', lambda, 'g', g, 'alpha', alpha);

num_iter = 30;
vfe = zeros(num_iter,1); 
% grad = zeros(num_iter,1);
% lambda = 10.^(linspace(-11,-9, 100));
% beta = ones(L_h_true, 1);
% alpha = ones(L_h_true, 1);
% ind = linspace(1,100,200);
% j = 1;
% for i = ind
%     beta(1) = i;
% %     vfe(j) = algo.computeVFE('beta', beta);
%     temp = algo.gradBetaTest(beta);
%     grad(j) = temp(1);
%     j = j + 1;
% end
% 
% plot(ind,grad)
vfe(1) = algo.computeVFE();
vfe(1)
t = 1:L_h;
algo.estimateH();
plot(t,algo.h_hat/max(abs(algo.h_hat)),t,h/max(h), 1:L_h_true, algo.alpha/max(algo.alpha),  1:L_h_true, algo.beta/max(algo.beta))
pause

tic
for i = 2:num_iter
%     algo.updateBeta('dichotomie', 1);
%     algo.updateAlpha('dichotomie', 1);
%     algo.updateBeta('dichotomie', 1);
    algo.updateG('conv');
    if mod(i,3) == 0
%         algo.updateLambda('dichotomie', 1);
%         algo.updateG('conv');
    end
% %     algo.lambda
    
    vfe(i) = algo.computeVFE();
    vfe(i)
    algo.estimateH();

%     plot(t,algo.h_hat/max(abs(algo.h_hat)),t,h/max(h), 1:L_h_true, algo.alpha/max(algo.alpha),  1:L_h_true, algo.beta/max(algo.beta))
%     pause;
end
toc
algo.estimateH();
plot(t,algo.h_hat/max(abs(algo.h_hat)),t,h/max(h), 1:L_h_true, algo.alpha/max(algo.alpha),  1:L_h_true, algo.beta/max(algo.beta))
