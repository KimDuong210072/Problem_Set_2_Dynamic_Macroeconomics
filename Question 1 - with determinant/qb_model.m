%% File Info.

%{

    model.m
    -------
    This code sets up the model.

%}

%% Model class.

classdef qb_model
    methods(Static)
        %% Set up structure array for model parameters and set the simulation parameters.
        function par = setup()            
            %% Structure array for model parameters.
            par = struct();
            
            %% Load Data
            % Load the data
            muc4a = readtable('muc4a.csv');
            muc123a = readtable('muc123a.csv');
            hh_expe = readtable('hhexpe08.csv'); 
            gold_data = readtable('muc5b4.csv'); 
            
            % Merge tables on keys: tinh, diaban, hoso, matv
            merged = outerjoin(muc4a, muc123a, ...
                'Keys', {'tinh','huyen','xa', 'diaban', 'hoso', 'matv'}, ...
                'MergeKeys', true);
            merged = outerjoin(merged, hh_expe, ...
                    'Keys', {'tinh','huyen','xa', 'diaban', 'hoso'}, ...
                    'MergeKeys', true);
            % Merge with gold purchase data
            merged = outerjoin(merged, gold_data, ...
                'Keys', {'tinh','huyen','xa','diaban','hoso'}, ...
                'MergeKeys', true);
            % Keep only household heads who are male
            is_male_head = merged.m1ac3 == 1 & merged.m1ac2 == 1;
         
            % Keep only valid income entries
            valid_income = merged.m4ac11 > 0 & ~ismissing(merged.m4ac11);
            
            % Filter the merged table
            filtered = merged(is_male_head & valid_income, :);
            
            % Compute log income
            i = filtered.m4ac11;
            par.i = i ;
            log_income = log(filtered.m4ac11);
            par.income = log_income;
            % Group by age and compute mean log income
            ages = filtered.m1ac5;
            [G, age_values] = findgroups(ages);
            mean_log_income = splitapply(@mean, log_income, G);
            
            % Exponentiate to get Gt
            Gt = exp(mean_log_income);
            
            % Display results
            results = table(age_values, mean_log_income, Gt, ...
                'VariableNames', {'Age', 'MeanLogIncome', 'Gt'});
            disp(results);
            
            % Optional: save to CSV
            writetable(results, 'Gt_by_age.csv');

            %% Load data for wealth and consumption
            par.c = filtered .riceexp + filtered .educex_2 + filtered .hlthex_1 + filtered .waterexp + filtered .elecexp;
            
            par.gold = filtered.m5b4c2_3; % Store gold purchases
            par.gold(isnan(par.gold)) = 0; % Replace missing with 0

            par.a0 = par.income - par.c - par.gold;
             %% Simulation parameters.
            par.seed = 2025; % Seed for simulation
            par.T = 60; %Last Period
            par.TT = 61; % Number of time periods.
            par.t_r = 41; % Retirement period
            par.NN = min(3000, length(par.a0)); 
            %% Preferences.
            par.beta = 0.94; % Discount factor
            %beta_values = [0.90, 0.92, 0.94, 0.96];
           
            % Store results
            %consumption_profiles = zeros(length(beta_values), par.T);
            %wealth_profiles = zeros(length(beta_values), par.T);
            
            % Loop over beta values
            %for i = 1:length(beta_values)
                % Setup and configure parameters
                %par = q_model.setup();
                %par.beta = beta_values(i);
            %end
            par.sigma = 2.00; % CRRA
            par.y_bar = 3.50; % Labor income
            par.k = 0.3;
            par.r = 0.15;
            par.phi = 0.005; % Annual appreciation of gold
            par.sigma_eps = 0.07; % Std. dev of productivity shocks
            par.rho = 0.85; % Persistence of AR(1) process
            par.mu = 0.0; % Intercept of AR(1) process
            
            %% Store G_t Data in Model Parameters
            par.Gt_table = results; % Save G_t_data in the parameter structure
            
            %% Handle Age Groups
            par.age_groups = unique(results.Age); % Extract unique age groups
            par.age_index = NaN(par.T, 1);

            % Loop through each model period and map to actual age index
            for t = 1:par.T
                [~, age_match] = min(abs(par.age_groups - t));  % Find closest age group
                par.age_index(t) = age_match;
            end
        
            %% Assertions (Error Checks)
            assert(par.beta > 0 && par.beta < 1.00,'Discount factor should be between 0 and 1.\n')
            assert(par.sigma >= 1,'CRRA should be at least 1.\n')
            assert(par.y_bar > 0,'Labor income should be larger than 0\n')
            assert(par.k > 0 && par.k < 1.00,'Pension fraction should be between 0 and 1.\n')
            assert(par.sigma_eps > 0,'The standard deviation of the shock must be positive.\n')
            assert(abs(par.rho) < 1,'The persistence must be less than 1 in absolute value so that the series is stationary.\n')

            

        end
        
        %% Generate state grids.
        
        function par = gen_grids(par)
            %% Cake grid.
             
            par.alen = 30; % Grid size for a.
            par.amax = 15.00; % Upper bound for a.
            par.amin = 0.00; % Minimum a.
            
            assert(par.alen > 5,'Grid size for a should be positive and greater than 5.\n')
            assert(par.amax > par.amin,'Minimum a should be less than maximum value.\n')
            
            par.agrid = linspace(par.amin,par.amax,par.alen)'; % Equally spaced, linear grid for a (and a').
                       
            %% Discretized productivity process.
            par.ylen = 30; % Grid size for yt.      
            par.m = 3; % Scaling parameter for Tauchen.
            assert(par.m > 0,'Scaling parameter for Tauchen should be positive.\n')
            assert(par.ylen > 3,'Grid size for y should be positive and greater than 3.\n')

            [ygrid,pmat] = qb_model.tauchen(par.mu,par.rho,par.sigma_eps,par.ylen,par.m); % Tauchen's Method to discretize the AR(1) process for log productivity.
            par.ygrid = exp(ygrid); % The AR(1) is in logs so exponentiate it to get y.
            par.pmat = pmat; % Transition matrix.

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
        
       %% Utility function.
        
        function u = utility(c,par)
            %% CRRA utility
            
            if par.sigma == 1
                u = log(c); % Log utility.
            else
                u = (c.^(1-par.sigma))./(1-par.sigma); % CRRA utility.
            end
                        
        end
        
    end
end