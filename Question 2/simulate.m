%% File Info.
%{
    simulate.m
    ----------
    This code simulates the model for multiple firms, including productivity
    and price shocks.
%}

%% Simulate class.
classdef simulate
    methods(Static)
        function sim = firm_dynamics(par, sol)
            %% Set up
            kgrid = par.kgrid;
            Agrid = par.Agrid;
            pgrid = par.pgrid;

            vpol = sol.v;
            kpol = sol.k;
            ipol = sol.i;
            rpol = sol.r;
            epol = sol.e;
            ppol = sol.p;

            T = par.T;
            N = par.N;
            burn = T; % Burn-in

            %% Containers for simulation results
            Asim = zeros(2*T, N);
            Psim = zeros(2*T, N);
            vsim = zeros(2*T, N);
            ksim = zeros(2*T, N);
            isim = zeros(2*T, N);
            rsim = zeros(2*T, N);
            esim = zeros(2*T, N);
            psim = zeros(2*T, N);

            %% Initialization
            rng(par.seed);  % Set random seed for reproducibility

            % Stationary distributions
            pmat0 = par.pmat^1000;
            pmat0 = pmat0(1,:);
            ppmat0 = par.pmat_p^1000;
            ppmat0 = ppmat0(1,:);

            % CDF matrices
            Acdf = cumsum(par.pmat, 2);   % For productivity transitions
            Pcdf = cumsum(par.pmat_p, 2); % For price transitions

            % Initial state indices
            k0_ind = randsample(par.klen, N, true);
            A0_ind = randsample(par.Alen, N, true, pmat0);
            P0_ind = randsample(par.plen, N, true, ppmat0);

            %% Initial values
            for i = 1:N
                Asim(1,i) = Agrid(A0_ind(i));
                Psim(1,i) = pgrid(P0_ind(i));
                vsim(1,i) = vpol(k0_ind(i), A0_ind(i), P0_ind(i));
                ksim(1,i) = kpol(k0_ind(i), A0_ind(i), P0_ind(i));
                isim(1,i) = ipol(k0_ind(i), A0_ind(i), P0_ind(i));
                rsim(1,i) = rpol(k0_ind(i), A0_ind(i), P0_ind(i));
                esim(1,i) = epol(k0_ind(i), A0_ind(i), P0_ind(i));
                psim(1,i) = ppol(k0_ind(i), A0_ind(i), P0_ind(i));

                % Draw next productivity and price indices
                uA = rand;
                A0_ind(i) = find(uA <= Acdf(A0_ind(i), :), 1);

                uP = rand;
                P0_ind(i) = find(uP <= Pcdf(P0_ind(i), :), 1);
            end

            %% Simulation loop
            for t = 2:(2*T)
                for i = 1:N
                    % Find nearest k index
                    [~, kt_ind] = min(abs(kgrid - ksim(t-1, i)));

                    % Record simulated states and policies
                    Asim(t,i) = Agrid(A0_ind(i));
                    Psim(t,i) = pgrid(P0_ind(i));
                    vsim(t,i) = vpol(kt_ind, A0_ind(i), P0_ind(i));
                    ksim(t,i) = kpol(kt_ind, A0_ind(i), P0_ind(i));
                    isim(t,i) = ipol(kt_ind, A0_ind(i), P0_ind(i));
                    rsim(t,i) = rpol(kt_ind, A0_ind(i), P0_ind(i));
                    esim(t,i) = epol(kt_ind, A0_ind(i), P0_ind(i));
                    psim(t,i) = ppol(kt_ind, A0_ind(i), P0_ind(i));

                    % Draw next productivity index
                    uA = rand;
                    A0_ind(i) = find(uA <= Acdf(A0_ind(i), :), 1);

                    % Draw next price index
                    uP = rand;
                    P0_ind(i) = find(uP <= Pcdf(P0_ind(i), :), 1);
                end
            end

            %% Collect results after burn-in
            sim = struct();
            sim.Asim = Asim(burn+1:end, :);
            sim.Psim = Psim(burn+1:end, :);
            sim.vsim = vsim(burn+1:end, :);
            sim.ksim = ksim(burn+1:end, :);
            sim.isim = isim(burn+1:end, :);
            sim.rsim = rsim(burn+1:end, :);
            sim.esim = esim(burn+1:end, :);
            sim.psim = psim(burn+1:end, :);
        end
    end
end
