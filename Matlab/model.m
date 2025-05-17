%% File Info.

%{

    model.m
    -------
    This code sets up the model.

%}

%% Model class.

classdef model
    methods(Static)
        %% Set up structure array for model parameters and set the simulation parameters.
        
        function par = setup()            
            %% Structure array for model parameters.
            
            par = struct();
            
            %% Technology.

            par.beta = 0.96; % Discount factor.
            par.alpha = 0.6; % Capital's share of income.
            par.delta = 0.6; % Depreciation rate.
            
            assert(par.delta >= 0.0 && par.delta <= 1.0,'The depreciation rate should be from 0 to 1.\n')
            assert(par.beta > 0.0 && par.beta < 1.0,'Discount factor should be between 0 and 1.\n')
            assert(par.alpha > 0.0 && par.alpha < 1.0,'Capital share of income should be between 0 and 1. \n')

            %% Prices, Income, and Costs
            %par.p = 1.00; % Price of investment.
            par.gamma = 1.00; % Speed of adjustment; cost function coefficient.

            par.sigma_eps = 0.07; % Std. dev of productivity shocks.
            par.rho = 0.85; % Persistence of AR(1) process.
            par.mu = 0.0; % Intercept of AR(1) process.

            %par.sigma_epsp = 0.05;
            %par.rhop = 0.85;
            %par.mup = 0.0;
            
            par.w = 5.0; 
            par.x = 2.0; 

            assert(par.gamma >= 0.0,'The cost function coefficient should be non-negative.\n')
            %assert(par.p > 0.0,'The price of investment should be positive.\n')

            assert(par.sigma_eps > 0,'The standard deviation of the shock must be positive.\n')
            assert(abs(par.rho) < 1,'The persistence must be less than 1 in absolute value so that the series is stationary.\n')
            %assert(par.sigma_epsp > 0,'The standard deviation of the shock must be positive.\n')
            %assert(abs(par.rhop) < 1,'The persistence must be less than 1 in absolute value so that the series is stationary.\n')

            %% Simulation parameters.

            par.seed = 2025; % Seed for simulation.
            par.T = 300; % Number of time periods.

        end
        
        %% Generate state grids.
        
        function par = gen_grids(par)
            %% Capital grid.

            par.klen = 300; % Grid size for a.
            par.kmax = 30.0; % Upper bound for a.
            par.kmin = 1e-4; % Minimum a.
            
            assert(par.klen > 5,'Grid size for k should be positive and greater than 5.\n')
            assert(par.kmax > par.kmin,'Minimum k should be less than maximum value.\n')
            
            par.kgrid = linspace(par.kmin,par.kmax,par.klen)'; % Equally spaced, linear grid for a and a'.
                
            %% Discretized income process.
                  
            par.Alen = 7; % Grid size for y.
            par.m = 3; % Scaling parameter for Tauchen.
            
            assert(par.Alen > 3,'Grid size for A should be positive and greater than 3.\n')
            assert(par.m > 0,'Scaling parameter for Tauchen should be positive.\n')
            
            [Agrid,pmat] = model.tauchen(par.mu,par.rho,par.sigma_eps,par.Alen,par.m); % Tauchen's Method to discretize the AR(1) process for log productivity.
            par.Agrid = exp(Agrid); % The AR(1) is in logs so exponentiate it to get A.
            par.pmat = pmat; % Transition matrix.

            par.plen = 8;
            assert(par.plen > 3,'Grid size for A should be positive and greater than 3.\n')
            %assert(par.m > 0,'Scaling parameter for Tauchen should be positive.\n')
            [pgrid, pmatP] = model.tauchen(par.mu, par.rho, par.sigma_eps, par.plen, par.m);
            par.pgrid = exp(pgrid);
            par.pmat_p = pmatP;
        
        end
        
        %% Tauchen's Method
        
        function [y,pi] = tauchen(mu,rho,sigma,N,m)
            %% Construct equally spaced grid.
        
            ar_mean = mu/(1-rho); % The mean of a stationary AR(1) process is mu/(1-rho).
            ar_sd = sigma/((1-rho^2)^(1/2)); % The std. dev of a stationary AR(1) process is sigma/sqrt(1-rho^2)
            
            y1 = ar_mean-(m*ar_sd); % Smallest grid point is the mean of the AR(1) process minus m*std.dev of AR(1) process.
            yn = ar_mean+(m*ar_sd); % Largest grid point is the mean of the AR(1) process plus m*std.dev of AR(1) process.
            
	        y = linspace(y1,yn,N); % Equally spaced grid.
            d = y(2)-y(1); % Step size.
	        
	        %% Compute transition probability matrix from state j (row) to k (column).
        
            ymatk = repmat(y,N,1); % States next period.
            ymatj = mu+rho*ymatk'; % States this period.
        
	        pi = normcdf(ymatk,ymatj-(d/2),sigma) - normcdf(ymatk,ymatj+(d/2),sigma); % Transition probabilities to state 2, ..., N-1.
	        pi(:,1) = normcdf(y(1),mu+rho*y-(d/2),sigma); % Transition probabilities to state 1.
	        pi(:,N) = 1 - normcdf(y(N),mu+rho*y+(d/2),sigma); % Transition probabilities to state N.
	        
        end
        
        
        %% Revenue function.
        
        function output = production(A,k,par)
            %% Revenue function.
            
            output = A.*k.^par.alpha; % Cobb-Douglas production.
                        
        end
        
        %% Cost function.
        
        function cost = total_cost(k,p,par)
            %% Convex adjustment cost.
            
            invest = par.kgrid-(1-par.delta).*k;
            adj_cost = (par.gamma/2).*((invest./k).^2).*k; % Convex adjustment cost.
            cost = adj_cost + p*invest; % Total investment cost.
                        
        end
    end
end