function [h,g] = synthesiseH(L_h, L_g, abs, lambda, var_h)
    
%     abs = 1e-4;
    L_h_true = L_h - L_g + 1;
    e_u = exp(-abs*((1:L_h)))';
    v = (1:L_h_true)';
    PI = gamrnd(lambda*v.^2,1./v);
%     g = conv(100*randn(L_g-2,1).*exp(-(1:L_g-2)/50)', [-1,2,-1]);
%     g = conv([1;-0.5; zeros(L_g-4,1)], [-1,2,-1]);
    g = [1; zeros(L_g-1,1)];
    b = conv(g, PI, 'full');
    w = randn(L_h, 1)*var_h;
    h = e_u.*b + w;

end