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
abs = 3*log(10)/(rt60*fe);
var_h = 1e-8;
L_g = 200;
L_h = 1e4;
L_h_true = L_h - L_g + 1;
beta = (1:L_h_true)';
alpha = (lambda*(1:L_h_true).^2)';
[h,g] = synthesiseH(L_h, L_g, abs, lambda, var_h);

algo = VEM('h', h, 'fe', fe, 'var_h', var_h, 'L_g', L_g, 'lambda', lambda, 'abs', abs, 'alpha', alpha);%, 'lambda', lambda, 'g', g, 'alpha', alpha);

num_iter = 2;
vfe = zeros(num_iter,1); 
% lambda = 10.^(linspace(-11,-9, 100));
% beta = ones(L_h_true, 1);
% alpha = ones(L_h_true, 1);
% ind = linspace(0.001,1000,1000);
% j = 1;
% for i = ind
%     alpha(100) = i;
%     vfe(j) = algo.computeVFE('alpha', alpha);
%     j = j + 1;
% end
% 
% semilogy(vfe)
vfe(1) = algo.computeVFE();
vfe(1)
t = 1:L_h;
tic
for i = 2:num_iter
%     algo.updateAlpha('dichotomie', 1);
    algo.updateBeta('dichotomie', 1);
    if mod(i,3) == 0
        algo.updateLambda('dichotomie', 1);
%         algo.updateG('conv');
    end
% %     algo.lambda
    
    vfe(i) = algo.computeVFE();
    vfe(i)
    algo.estimateH();

%     plot(t,algo.h_hat/max(algo.h_hat),t,h/max(h), 1:L_h_true, algo.alpha/max(algo.alpha),  1:L_h_true, algo.beta/max(algo.beta))
%     pause;
end
toc
algo.estimateH();
plot(t,algo.h_hat/max(algo.h_hat),t,h/max(h), 1:L_h_true, algo.alpha/max(algo.alpha),  1:L_h_true, algo.beta/max(algo.beta))
