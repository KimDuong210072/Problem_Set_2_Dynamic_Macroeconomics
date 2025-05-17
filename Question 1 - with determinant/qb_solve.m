%% Solve class.

classdef qb_solve
    methods(Static)
        %% Solve the model using BI. 
        function sol = cs_qb_model_fin(par)  
            sol = struct();
            
            %% Extract parameters
            T = par.T;
            tr = par.t_r;
            beta = par.beta;
            r = par.r;
            kappa = par.k;
        
            alen = par.alen;
            agrid = par.agrid;
        
            ylen = par.ylen;
            ygrid = par.ygrid;
            pmat = par.pmat;
        
            %age_groups = par.age_groups;
            %age_index = par.age_index;
            Gt_log = log(par.Gt_table.Gt);  % take log so we can exponentiate later
            %rho = par.rho;
            
            %% Initialize containers
            v1 = nan(alen, T, ylen);
            a1 = nan(alen, T, ylen);
            c1 = nan(alen, T, ylen);
            g1 = nan(alen, glen, T, ylen);

            amat = repmat(agrid, 1, ylen);
            ymat = repmat(ygrid, alen, 1);
        
            fprintf('------------Solving from the Last Period of Life.------------\n\n')
        
            for age = 1:T
                t = T - age + 1;  % Actual time index
        
                if t == T
                    % Last period: consume everything
                    c1(:, T, :) = amat + kappa * ymat;
                    a1(:, T, :) = 0;
                    v1(:, T, :) = qb_model.utility(c1(:, T, :), par);
                else
                    for i = 1:ylen
                        % Determine income y_t based on age
                        if t >= tr
                            % Retired: constant income
                            yt = kappa * ygrid(i);
                            ev = v1(:,T-age+2,i);
                        else
                            % Working: stochastic income
                            age_idx = par.age_index(t);
                            yt = Gt_log(age_idx)*ygrid(i);  % log income process
                            ev = squeeze(v1(:,T-age+2,:))*pmat(i,:)';
                        end
        
                        for p = 1:alen  % current asset level (a)
                        val = -inf(alen, 1);  % Initialize value function for each a'
                        ct_all = agrid(p) + yt - (agrid / (1 + r));  % All possible future asset levels
                    
                        for aprime = 1:alen
                            ct = ct_all(aprime);
                            if ct > 0
                                u = qb_model.utility(ct, par);
                                val(aprime) = u + beta * ev(aprime);
                            end
                        end
                            [vmax, ind] = max(val);
        
                            v1(p, t, i) = vmax;
                            c1(p, t, i) = ct_all(ind);
                            a1(p, t, i) = agrid(ind);
                        end
                    end
                end
        
                if mod(t, 5) == 0
                    fprintf('Age: %d.\n', t);
                end
            end
        
            fprintf('------------Life Cycle Problem Solved.------------\n')
        
            sol.c = c1;
            sol.a = a1;
            sol.v = v1;
        end
    end
end    