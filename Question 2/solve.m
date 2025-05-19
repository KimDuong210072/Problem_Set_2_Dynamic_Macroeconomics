%% File Info.
%{
    solve.m
    -------
    This code solves the model using VFI for big and small firms.
%}

%% Solve class.

classdef solve
    methods(Static)
        %% Solve the model using VFI.
        
        function sol = firm_problem(par)
            %% Structure array for model solution.
            
            sol = struct();
            
            %% Model parameters, grids and functions.
            
            beta = par.beta;
            delta = par.delta;
            
            klen = par.klen;
            kgrid = par.kgrid;
            
            Alen = par.Alen;
            Agrid = par.Agrid;
            pmat = par.pmat;
            
            plen = par.plen;
            pgrid = par.pgrid;
            pmatp = par.pmat_p;
            
            w = par.w;
            x = par.x;
            
            %% Initialize value and policy functions.
            V = zeros(klen, Alen, plen);
            Vnew = zeros(klen, Alen, plen);
            policy_kprime = zeros(klen, Alen, plen);
            policy_invest = zeros(klen, Alen, plen);
            revenue_fun = zeros(klen, Alen, plen);
            profit_fun = zeros(klen, Alen, plen);
            cost_fun = zeros(klen, Alen, plen);
            
            crit = 1e-6;
            maxiter = 10000;
            diff = 1;
            iter = 0;

            if isfield(par, 'data')
                firm_type = unique(par.data.group);
                fprintf('Solving for firm group: %s\n', firm_type);
            else
                fprintf('Solving for firm group: UNKNOWN\n');
            end

            fprintf('------------Beginning Value Function Iteration.------------\n');

            while diff > crit && iter < maxiter
                for iA = 1:Alen
                    A = Agrid(iA);
                    for iP = 1:plen
                        p = pgrid(iP);
                        for ik = 1:klen
                            k0 = kgrid(ik);

                            % Value for each k1 choice
                            Vchoices = zeros(klen, 1);
                            for ikp = 1:klen
                                kp = kgrid(ikp);

                                % Revenue and cost
                                rev = model.production(A, k0, par);
                                cst = model.total_cost(k0, kp, p, par);

                                % Expected future value
                                EV = 0;
                                for iA2 = 1:Alen
                                    for iP2 = 1:plen
                                        prob = pmat(iA, iA2) * pmatp(iP, iP2);
                                        EV = EV + prob * V(ikp, iA2, iP2);
                                    end
                                end

                                Vchoices(ikp) = rev - cst + beta * EV;
                            end

                            % Maximize over future capital choices
                            [Vnew(ik, iA, iP), best_idx] = max(Vchoices);
                            best_kp = kgrid(best_idx);
                            policy_kprime(ik, iA, iP) = best_kp;

                            % Record policy functions
                            investment = best_kp - (1 - delta) * k0;
                            revenue_fun(ik, iA, iP) = model.production(A, k0, par);
                            cost_fun(ik, iA, iP) = model.total_cost(k0, best_kp, p, par);
                            policy_invest(ik, iA, iP) = investment;
                            profit_fun(ik, iA, iP) = revenue_fun(ik, iA, iP) - cost_fun(ik, iA, iP);
                        end
                    end
                end

                diff = max(abs(Vnew(:) - V(:)));
                V = Vnew;
                iter = iter + 1;
                fprintf('Iteration %d, diff = %.8f\n', iter, diff);
            end

            fprintf('Converged in %d iterations.\n', iter);
            fprintf('------------End of Value Function Iteration.------------\n');

            %% Store solution
            sol.v = V;
            sol.k = policy_kprime;
            sol.i = policy_invest;
            sol.r = revenue_fun;
            sol.e = cost_fun;
            sol.p = profit_fun;
        end
    end
end
