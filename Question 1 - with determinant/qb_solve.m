classdef qb_solve
    methods(Static)
        function sol = cs_qb_model_fin(par)
            sol = struct();

            %% Extract parameters
            T = par.T;
            tr = par.t_r;
            beta = par.beta;
            r = par.r;
            kappa = par.k;
            phi = par.phi;

            alen = par.alen;
            agrid = par.agrid;

            ylen = par.ylen;
            ygrid = par.ygrid;
            pmat = par.pmat;

            Gt = par.Gt_table.Gt;  % Correct: use level value
            gold0 = mean(par.gold);  % Assume initial gold value

            %% Initialize containers
            v1 = nan(alen, T, ylen);
            a1 = nan(alen, T, ylen);
            c1 = nan(alen, T, ylen);

            amat = repmat(agrid, 1, ylen);
            ymat = repmat(ygrid, alen, 1);

            fprintf('------------Solving from the Last Period of Life.------------\n\n')

            for age = 1:T
                t = T - age + 1;  % Backward time index

                % Gold evolution
                if t == 1
                    gold_value = -gold0;  % Purchase gold
                else
                    gold_value = gold0 * (1 + phi)^(t - 1);  % Accumulated value
                end

                if t == T
                    % Terminal period: consume all available resources
                    c1(:, T, :) = amat + kappa * ymat + gold_value;
                    a1(:, T, :) = 0;
                    v1(:, T, :) = qb_model.utility(c1(:, T, :), par);
                else
                    for i = 1:ylen
                        % Income determination
                        if t >= tr
                            yt = kappa * ygrid(i);  % Pension income
                            ev = v1(:, t + 1, i);    % No income shocks
                        else
                            age_idx = par.age_index(t);
                            yt = Gt(age_idx) * ygrid(i);  % Level income
                            ev = squeeze(v1(:, t + 1, :)) * pmat(i, :)';  % Expectation
                        end

                        for p = 1:alen
                            val = -inf(alen, 1);
                            available = (1 + r) * (agrid(p) + yt + gold_value);  % Budget

                            ct_all = available - agrid;
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
