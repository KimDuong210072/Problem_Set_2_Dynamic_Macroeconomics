%% File Info.

%{

    solve.m
    -------
    This code solves the model.

%}

%% Solve class.

classdef solve
    methods(Static)
        %% Solve the model using VFI. 
        
        function sol = firm_problem(par)            
            %% Structure array for model solution.
            
            sol = struct();
            
            %% Model parameters, grids and functions.
            
            beta = par.beta; % Discount factor.
            delta = par.delta; % Depreciation rate.

            klen = par.klen; % Grid size for k.
            kgrid = par.kgrid; % Grid for k (state and choice).

            Alen = par.Alen; % Grid size for A.
            Agrid = par.Agrid; % Grid for A.
            pmat = par.pmat; % Transition matrix for A.

            plen = par.plen; % Length of price grid
            pgrid = par.pgrid; % Grid for price
            pmatp = par.pmat_p; % Transition matrix for Price.

            w = par.w;
            x = par.x;

            %% Value Function iteration.

            v0 = zeros(klen,Alen,plen); % Guess of value function is zero profit.

            v1 = nan(klen,Alen,plen); % Container for V.
            k1 = nan(klen,Alen,plen); % Container for K'.
            i1 = nan(klen,Alen,plen); % Container for i.
            r1 = nan(klen,Alen,plen); % Container for revenue.
            e1 = nan(klen,Alen,plen); % Container for investment expenditure.
            p1 = nan(klen,Alen,plen); % Container for profit.

            crit = 1e-6;
            maxiter = 10000;
            diff = 1;
            iter = 0;
            
            fprintf('------------Beginning Value Function Iteration.------------\n\n')
            
            while diff > crit && iter < maxiter % Iterate on the Bellman Equation until convergence.
                
                for p = 1:klen % Loop over the K-states.
                    for j = 1:Alen % Loop over the A-states.
                        for h = 1:plen
                            % Macro variables.
                            rev = model.production(Agrid(j), kgrid(p), par);  % Current output
                            vall = -inf(par.klen, 1);         % Value function over choices
                            invest_vec = nan(par.klen, 1);    % Investment per k′
                            prof_vec = nan(par.klen, 1);      % Profit per k′
                            
                            for i = 1:par.klen
                                k_prime = kgrid(i);
                                invest = k_prime - (1 - par.delta) * kgrid(p);
                            
                                if invest <= 0
                                    continue
                                end
                            
                                adj_cost = (par.gamma / 2) * ((invest / kgrid(p))^2) * kgrid(p);
                                inv_cost = adj_cost + pgrid(h) * invest;
                            
                                profit = rev - par.x * par.w - inv_cost;
                            
                                ev = 0;
                                for jj = 1:par.Alen
                                    for hh = 1:par.plen
                                        ev = ev + v0(i, jj, hh) * pmat(j, jj) * pmatp(h, hh);
                                    end
                                end
                            
                                vall(i) = profit + par.beta * ev;
                                invest_vec(i) = invest;
                                prof_vec(i) = profit;
                            end
                            
                            [vmax, ind] = max(vall);
                            
                            v1(p, j, h) = vmax;
                            k1(p, j, h) = kgrid(ind);
                            i1(p, j, h) = invest_vec(ind);
                            r1(p, j, h) = rev;
                            e1(p, j, h) = inv_cost;  
                            p1(p, j, h) = prof_vec(ind); 

                        end
                    end
                end
                
                diff = norm(v1-v0,'fro'); % Check for convergence.
                v0 = v1; % Update guess of v.
                
                iter = iter + 1; % Update counter.
                
                % Print counter.
                if mod(iter,25) == 0
                    fprintf('Iteration: %d.\n',iter)
                end

            end
                
            fprintf('\nConverged in %d iterations.\n\n',iter)
            
            fprintf('------------End of Value Function Iteration.------------\n')
            
            %% Macro variables, value, and policy functions.
            
            sol.v = v1; % Firm value.
            sol.k = k1; % Capital policy function.
            sol.i = i1; % Investment policy function.
            sol.r = r1; % Revenue function.
            sol.e = e1; % Investment expenditure function.
            sol.p = p1; % Profit function.
            
        end
        
    end
end