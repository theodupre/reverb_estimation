classdef VEM < handle

    properties
        
        h; % observed variable : RIR
        h_hat % estimate of h
        fe; % sample rate
        L_h; % length of h
        
        %%%%%%%%%%%
        % Model parameters
        %%%%%%%%%%%
        g; % Microphone/Source IR
        gg; % test
        var_h; % variance of noise of h
        lambda; % mean of Poisson process
        L_g; % length of g
        L_h_true; % L_h - L_g + 1
        a; % absorption coefficient
        
        
        %%%%%%%%%%%
        % Variational parameters
        %%%%%%%%%%%
        alpha; 
        beta;
        
        %%%%%%%%%%
        % Others
        %%%%%%%%%%
        u; % Time variable
        v; % Time variable
        PI; % distribution to be approximated
        b; % g * PI
        e_u; % exp(-au)/var_h
        e_2u; % exp(-2au)/var_h
        g_p;
        beta_p;
        
    end
    
    methods
        
        function this = VEM(varargin)
            this.Set(varargin{:});
        end
        
        function Set(this,varargin)
            st_input = struct(varargin{:});
            
            if isfield(st_input, 'h')
                this.h = st_input.h;
            end
            
            if isfield(st_input, 'var_h')
                this.var_h = st_input.var_h;
            end
            
            if isfield(st_input, 'fe')
                this.fe = st_input.fe;
            end
            
            if isfield(st_input, 'L_g')
                this.L_g = st_input.L_g;
            end
            
            this.L_h = length(this.h);
            this.L_h_true = this.L_h - this. L_g + 1;
            this.u = (1:this.L_h)';
            this.v = (1:this.L_h_true)' ;
            
            if isfield(st_input, 'alpha')
                this.alpha = st_input.alpha;
            else
                this.alpha = ones(this.L_h_true,1)*1000;
            end
            
            if isfield(st_input, 'beta')
                this.beta = st_input.beta;
            else
                this.beta = ones(this.L_h_true,1)*10;
            end
            %%%%%%
            % Initialization of model parameters
            %%%%%%
            
            if isfield(st_input, 'abs')
                this.a = st_input.abs;
            else
                h_log = log(abs(this.h + eps));
                reverb_ind = (floor(0.2*this.fe):floor(0.3*this.fe))';
                h_reverb = h_log(reverb_ind) + eps;
                fit = polyfit(reverb_ind, h_reverb, 1);
                this.a = -fit(1);
            end
            
            this.b = this.h.*exp(this.a*this.u);
            
            if isfield(st_input, 'lambda')
                this.lambda = st_input.lambda;
            else
                [~,this.lambda] = lpc(this.b, 50);
            end
            
            if isfield(st_input, 'g')
                this.g = st_input.g;
            else
                [A,~] = lpc(this.b, 50);
                this.g = filter(1,A,[1, zeros(1,this.L_g-1) ])';
                this.g = [1;zeros(this.L_g-1,1)];
            end

            this.PI = randn(this.L_h_true,1)*this.lambda.^2;
            this.e_u = exp(-this.a*this.u)/this.var_h;
            this.e_2u = exp(-2*this.a*this.u)/this.var_h;
            this.g_p = this.g.*exp(-this.a*this.u(1:this.L_g));
            this.beta_p = this.beta.*exp(this.a*this.v);
        
        end
        
        function vfe = computeVFE(this, param, x)
            if nargin > 1
                if strcmp(param, 'lambda')
                    pi_1 = x*this.v.^2.*log(this.v) - gammaln(x*this.v.^2);
                    pi_2 = x*this.v.^2.*(psi(0, this.alpha) - log(this.beta));

                    vfe = sum(pi_1 + pi_2);
                elseif strcmp(param, 'lambda_square')
                    pi_1 = x^2*this.v.^2.*log(this.v) - gammaln(x^2*this.v.^2);
                    pi_2 = x^2*this.v.^2.*(psi(0, this.alpha) - log(this.beta));

                    vfe = sum(pi_1 + pi_2);
                
                elseif strcmp(param, 'g')                
                    expect_term = (1./(this.e_u*this.var_h).*this.h - conv(x, this.alpha./this.beta, 'full')).^2;
                    var_term = conv(x, this.alpha./(this.beta.^2), 'full');   

                    vfe = -sum(this.e_2u/2.*(expect_term + var_term));

                elseif strcmp(param, 'beta')
                    x_p = x./(this.e_u(1:this.L_h_true).*this.var_h);
                    expect_term = sum((this.h - conv(this.g_p, this.alpha./x_p, 'full')).^2);
                    var_term = sum(conv(this.g_p.^2, this.alpha./(x_p.^2), 'full'));
                    pi_term = - log(x).*this.lambda.*this.v.^2  - this.alpha.*this.v./x;

                    vfe = - 1/(2*this.var_h)*(expect_term + var_term) + sum(pi_term);

                elseif strcmp(param, 'alpha')
                    expect_term = sum((this.h - conv(this.g_p, x./this.beta_p, 'full')).^2);
                    var_term = sum(conv(this.g_p.^2, x./(this.beta_p.^2), 'full'));
                    pi_term = gammaln(x) + (this.lambda.*this.v.^2 - x).*(psi(0,x) - log(this.beta)) - x.*(log(this.beta) + this.v./this.beta - 1);

                    vfe = - 1/(2*this.var_h)*(expect_term + var_term) + sum(pi_term);

                end
            else

                var_h_term = -this.L_h/2*log(2*pi*this.var_h);
                expect_term = (1./(this.e_u*this.var_h).*this.h - conv(this.g, this.alpha./this.beta, 'full')).^2;
                var_term = conv(this.g.^2, this.alpha./(this.beta.^2), 'full');
                
                norm_expect = sum((this.h - conv(this.g_p, this.alpha./this.beta_p,'full')).^2);
                norm_var = sum(this.g_p.^2).*sum(this.alpha./(this.beta_p.^2));
                
                likelihood_term = var_h_term - sum(this.e_2u/2.*(expect_term + var_term));
                likelihood_term2 = var_h_term - 1/(2*this.var_h)*(norm_expect + norm_var);
                
                pi_1 = this.lambda*this.v.^2.*log(this.v) + gammaln(this.alpha) - gammaln(this.lambda*this.v.^2);
                pi_2 = (this.lambda*this.v.^2 - this.alpha).*(psi(0,this.alpha) - log(this.beta));
                pi_3 = this.alpha.*(log(this.beta) + this.v./this.beta - 1);

                pi_term = sum(pi_1 + pi_2 - pi_3);

                vfe = likelihood_term + pi_term;
                vfe = likelihood_term2 + pi_term;
            end
        end
        
        function [] = updateLambda(this, mode, num_iter)
        
            if strcmp(mode, 'matlab')
                f = @(x) this.computeVFE('lambda', x);
                this.lambda = fminbnd(f, 0, 100);
            end
            
            if strcmp(mode, 'matlab_square')
                f = @(x) this.computeVFE('lambda_square', x);
                this.lambda = fminsearch(f, sqrt(this.lambda))^2;
            end
            
            function grad = gradLambda(param)
                grad = sum(this.v.^2.*(log(this.v) - psi(0, param*this.v.^2) + psi(0, this.alpha) - log(this.beta)));
            end
            
            function hess = hessLambda(param)
                hess = -sum(this.v.^4.*psi(1,param.*this.v.^2));
            end
            
            if strcmp(mode, 'dichotomie')
                down_limit = 1e-20;
                up_limit = 100;
                while (up_limit - down_limit)/up_limit > 1e-6
                   if  gradLambda((up_limit + down_limit)/2) > 0
                       down_limit = (up_limit + down_limit)/2;
                   else
                       up_limit = (up_limit + down_limit)/2;
                   end
                end
                this.lambda = (up_limit + down_limit)/2;
            
            end
            
        end

        function [] = updateG(this, mode, num_iter)
           
            if strcmp(mode, 'matlab')
                f = @(x) this.computeVFE('g', x);
                this.g = fminsearch(f,this.g);            
            end
            
            if strcmp(mode, 'newton')
                param = this.g;
                hess = -conv(this.e_2u, flip((this.alpha.^2 + this.alpha)./(this.beta.^2)), 'valid');
                mean(hess)
                for i = 1:num_iter
                    new_param = param - gradG(param)./hess;
                    param = new_param;
                end
                this.g = new_param;
            end
            
            if strcmp(mode, 'grad')
                toep_col = flip(this.alpha./this.beta);
                toep_col2 = flip(this.alpha./(this.beta.^2));
                diag_term = repmat(this.e_2u, 1, length(this.g)) ;
                
                toep_mat = toeplitz([toep_col; zeros(length(this.g)-1,1)], [1, ones(1,length(this.g) - 1)]);
                toep_mat2 = toeplitz([toep_col2; zeros(length(this.g)-1, 1)], [1, zeros(1,length(this.g)-1)]);
                
                M = toep_mat'*(diag_term.*toep_mat) + repmat(toep_mat2'*this.e_2u,1,length(this.g));
                C = toep_mat'*(this.e_u.*this.h);
                
                clear toep_col2 toep_mat diag_term
                
                this.g = C\M;
            end
            
            if strcmp(mode, 'conv')
               beta_p = this.beta./(this.e_u(1:this.L_h_true)*this.var_h);
               B = conv(this.h, flip(this.alpha)./flip(beta_p), 'full');
               B = B(1:this.L_g);
               c = sum(flip(this.alpha)./flip(beta_p.^2));
               gamma = conv(this.alpha./beta_p, flip(this.alpha)./beta_p, 'full');
               A = c*eye(this.L_g, this.L_g) + toeplitz(gamma(1:this.L_g));
               g_p = B\A;
               this.g = g_p'./(this.e_u(1:this.L_g)*this.var_h);
            end
        
            function grad = gradG(param)
                
                h_term = conv(this.e_u.*this.h, flip(this.alpha./this.beta), 'valid');
                conv_term = conv(this.e_2u.*conv(param, this.alpha./this.beta, 'full'), flip(this.alpha./this.beta), 'valid');
                var_term = param.*conv(this.e_2u, flip(this.alpha./(this.beta.^2)), 'valid');
                
                grad = h_term - conv_term - var_term;
                mean(grad)
            end
            
        end
        function updateAlpha(this, mode, num_iter)
            
            if strcmp(mode, 'newton')
                new_param = this.alpha;
                for i=1:num_iter
                    param = new_param;
                    new_param = param - gradAlpha(param)./hessAlpha(param);
                    if mod(i, 10) == 0
                       mean(new_param)
                    end
                end
                
                this.alpha = new_param;
            end
            
            if strcmp(mode, 'newton_exp')
                new_param = log(this.alpha);
                for i=1:num_iter
                    param = new_param;
                    new_param = param - gradLogAlpha(param)./hessLogAlpha(param);
                    if mod(i, 10) == 0
                       mean(exp(new_param))
                    end
                end
                
                this.alpha = exp(new_param);
                
            end
            
            if strcmp(mode, 'grad_desc')
                eta = 1e-4;
                new_param = sqrt(this.alpha);
                for i=1:num_iter
                    param = new_param;
                    new_param = param + eta*gradSqrtAlpha(param);
                end
                this.alpha = new_param.^2;
                
            end
            
            if strcmp(mode, 'newton_square')
                new_param = sqrt(this.alpha);
                for i=1:num_iter
                    param = new_param;
                    
                    if max(hessSqrtAlpha(param)) >= 0
                       new_param = param - gradSqrtAlpha(param)./(-abs(hessSqrtAlpha(param)));
                    else
                       new_param = param - gradSqrtAlpha(param)./hessSqrtAlpha(param);
                    end
                end
                
                this.alpha = new_param.^2;
                
            end
            
            function grad = gradAlpha(param)
                
                h_term = 1./this.beta.*conv(this.e_u.*this.h, flip(this.g), 'valid');
                expect_term = 1./this.beta.*conv(this.e_2u.*conv(this.g, param./this.beta, 'full'), flip(this.g), 'valid');
                var_term = 1./(this.beta.^2).*conv(this.e_2u/2, flip(this.g.^2), 'valid');
                pi_term = (this.lambda.*this.v.^2 - param).*psi(1, param) - this.v./this.beta + 1;
                
                grad = h_term - expect_term - var_term + pi_term;
                
            end
            
            function hess = hessAlpha(param)
               var_term = conv(this.e_2u, flip(this.g.^2), 'valid');
               pi_term = - psi(1, param) + (this.lambda.*this.v.^2 - param).*psi(2, param);
               
               hess = -1./(this.beta.^2).*var_term + pi_term;
            end
            
            function grad = gradLogAlpha(param)
                grad = gradAlpha(exp(param)).*exp(param);
            end
            
            function hess = hessLogAlpha(param)
                hess = hessAlpha(exp(param)).*exp(2*param) + gradAlpha(exp(param)).*exp(param);
            end
            
            function grad = gradSqrtAlpha(param)
                grad = 2*gradAlpha(param.^2).*param;
            end
            
            function hess = hessSqrtAlpha(param)
                hess = 4*hessAlpha(param.^2).*param.^2 + 2*gradAlpha(param.^2);
            end
        end
        
        function [] = updateBeta(this, mode, num_iter)
            
            if strcmp(mode, 'newton')
                new_param = this.beta;
                for i = 1:num_iter
                   param = new_param;
                   new_param = param - gradBeta(param)./hessBeta(param);
                end
                
                this.beta = new_param;
            end
            
            if strcmp(mode, 'newton_square')
                new_param = sqrt(this.beta);
                for i = 1:num_iter
                    param = new_param;
                    gradSqrtParam = gradSqrtBeta(param);
                    hessSqrtParam = hessSqrtBeta(param);
                    if max(hessSqrtParam)>=0 
                        new_param = param - gradSqrtParam./(-abs(hessSqrtParam));
                    else
                        new_param = param - gradSqrtParam./hessSqrtParam;   
                    end
                end
                
                this.beta = new_param.^2;
            end
            
            if strcmp(mode, 'grad_desc')
                eta = 1e-8;
                new_param = sqrt(this.beta);
                
                for i=1:num_iter
                    param = new_param;
                    new_param = param + eta*gradSqrtBeta(param);
                end
                
                this.beta = new_param.^2;
            
            end
            
            if strcmp(mode, 'newton_exp')
                new_param = log(this.beta);
                for i = 1:num_iter
                    param = new_param;
                    gradLogParam = gradBeta(exp(param)).*exp(param);
                    hessLogParam = hessBeta(exp(param)).*exp(2*param) + gradLogParam;
                    if max(hessLogParam)>=0 || min(hessLogParam)<=0 
                        warning('HessLogParam contains zero value, param not updated')
                    else
                        new_param = param - gradLogParam./hessLogParam;
                    end
                end
                
                this.beta = exp(new_param);
            end
            
            function grad = gradBeta(param)
                
                h_term = this.alpha./(param.^2).*conv(this.e_u.*this.h, flip(this.g), 'valid');
                expect_term = this.alpha./(param.^2).*conv(this.e_2u.*conv(this.g, this.alpha./param, 'full'), flip(this.g), 'valid');
                var_term =  this.alpha./(param.^3).*conv(this.e_2u, flip(this.g.^2), 'valid');
                pi_term = -this.lambda.*this.v.^2./param + this.v.*this.alpha./(param.^2);
                
                grad = -h_term + expect_term + var_term + pi_term;
                
            end
            
            function hess = hessBeta(param)
               
                h_term = 2*this.alpha./(param.^3).*conv(this.e_u.*this.h, flip(this.g), 'valid');
                expect_term = 2*this.alpha./(param.^3).*conv(this.e_2u.*conv(this.g, this.alpha./param, 'full'), flip(this.g), 'valid');
                quad_term = this.alpha.^2./(param.^4).*conv(this.e_2u, flip(this.g.^2), 'valid');
                var_term = 3*this.alpha./(param.^4).*conv(this.e_2u, flip(this.g.^2), 'valid');
                pi_term = this.lambda.*this.v.^2./(param.^2) - 2*this.alpha.*this.v./(param.^3);
                
                hess = h_term - expect_term - quad_term - var_term + pi_term;
               
            end
            
            function grad = gradSqrtBeta(param)
                grad = 2*gradBeta(param.^2).*param;
            end
            
            function hess = hessSqrtBeta(param)
                hess = 4*hessBeta(param.^2).*param.^2 + 2*gradBeta(param.^2);
            end
        end
            
        function [] = estimateH(this)
            
            this.PI = gamrnd(this.alpha, 1./this.beta);
            this.b = conv(this.PI, this.g, 'full'); 
            this.h_hat = this.e_u.*this.b + randn(this.L_h,1)*this.var_h;
        end
        
    end
    
end

