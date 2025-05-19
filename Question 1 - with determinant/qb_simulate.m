classdef qb_simulate
    methods(Static)
        function sim = lc(par, sol)
            % Set up
            agrid = par.agrid;
            apol = sol.a;
            cpol = sol.c;

            TT = par.TT;
            NN = par.NN;
            T = par.T;
            tr = par.t_r;

            kappa = par.k;
            ygrid = par.ygrid;
            pmat = par.pmat;

            phi = par.phi;
            gold0 = mean(par.gold);

            ysim = nan(TT, NN);
            Asim = nan(TT, NN);
            tsim = nan(TT, NN);
            csim = nan(TT, NN);
            usim = nan(TT, NN);
            goldsim = nan(TT, NN);

            rng(par.seed);

            pmat0 = pmat^100;
            cmat = cumsum(pmat, 2);

            y0_ind = randsample(par.ylen, NN, true, pmat0(1,:))';
            a0_ind = arrayfun(@(a0) find(abs(agrid - a0) == min(abs(agrid - a0)), 1), par.a0(:));
            t0_ind = randsample(T, NN, true)'; % Initial age
            yr = nan(NN, 1); % Retirement income

            for i = 1:NN
                age0 = t0_ind(i);
                tsim(1,i) = age0;

                if age0 == 1
                    goldsim(1,i) = -gold0;
                else
                    goldsim(1,i) = gold0 * (1 + phi)^(age0 - 1);
                end

                if age0 >= tr
                    yr(i) = ygrid(y0_ind(i));
                    ysim(1,i) = kappa * yr(i);
                else
                    ysim(1,i) = ygrid(y0_ind(i));
                end

                csim(1,i) = cpol(a0_ind(i), age0, y0_ind(i)) + goldsim(1,i);
                Asim(1,i) = apol(a0_ind(i), age0, y0_ind(i));
                usim(1,i) = qb_model.utility(csim(1,i), par);

                if age0 < tr - 1
                    draw = rand;
                    y0_ind(i) = find(draw <= cmat(y0_ind(i),:), 1, 'first');
                elseif age0 == tr - 1
                    yr(i) = ygrid(y0_ind(i));
                end
            end

            % Recursive simulation
            for j = 2:TT
                for i = 1:NN
                    age = tsim(j-1,i) + 1;

                    if age <= T
                        tsim(j,i) = age;
                        if age >= tr
                            yval = kappa * yr(i);
                        else
                            yval = ygrid(y0_ind(i));
                        end
                        ysim(j,i) = yval;

                        if age == 1
                            goldsim(j,i) = -gold0;
                        else
                            goldsim(j,i) = gold0 * (1 + phi)^(age - 1);
                        end

                        % Find closest index
                        [~, at_ind] = min(abs(agrid - Asim(j-1,i)));
                        if at_ind < 1 || at_ind > length(agrid)
                            at_ind = max(1, min(at_ind, length(agrid)));  % Clamp index
                        end

                        csim(j,i) = cpol(at_ind, age, y0_ind(i)) + goldsim(j,i);
                        Asim(j,i) = apol(at_ind, age, y0_ind(i));
                        usim(j,i) = qb_model.utility(csim(j,i), par);

                        if age < tr - 1
                            draw = rand;
                            y0_ind(i) = find(draw <= cmat(y0_ind(i),:), 1, 'first');
                        elseif age == tr - 1
                            yr(i) = ygrid(y0_ind(i));
                        end
                    end
                end
            end

            % Output struct
            sim = struct();
            sim.ysim = ysim;
            sim.Asim = Asim;
            sim.tsim = tsim;
            sim.csim = csim;
            sim.usim = usim;
            sim.goldsim = goldsim;
        end
    end
end
