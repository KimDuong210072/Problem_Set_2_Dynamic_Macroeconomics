%% File Info.

%{

    simulate.m
    ----------
    This code simulates the model.

%}

%% Simulate class.

classdef simulate
    methods(Static)
        %% Simulate the model. 
        
        function sim = firm_dynamics(par,sol)            
            %% Set up.
            
            kgrid = par.kgrid; % Capital today (state variable).
            Agrid = par.Agrid; % Productivity (state variable).

            vpol = sol.v; % Firm value.
            kpol = sol.k; % Policy function for capital.
            ipol = sol.i; % Policy function for investment.
            rpol = sol.r; % Optimal revenue.
            epol = sol.e; % Optimal total investment expenditure.
            ppol = sol.p; % Optimal profit.

            T = par.T; % Time periods.
            Asim = zeros(T*2,1); % Container for simulated productivity.
            vsim = zeros(T*2,1); % Container for simulated firm value.
            rsim = zeros(T*2,1); % Container for simulated output.
            ksim = zeros(T*2,1); % Container for simulated capital stock.
            isim = zeros(T*2,1); % Container for simulated investment.
            esim = zeros(T*2,1); % Container for simulated investment expenditure.
            psim = zeros(T*2,1); % Container for simulated profit.
            
            %% Begin simulation.
            
            rng(par.seed);

            pmat0 = par.pmat^1000;
            pmat0 = pmat0(1,:); % Stationary distribution.
            cmat = cumsum(par.pmat,2); % CDF matrix.

            k0_ind = randsample(par.klen,1); % Index for initial capital stock.
            A0_ind = randsample(par.Alen,1,true,pmat0); % Index for initial productivity.

            Asim(1) = Agrid(A0_ind); % Productivity in period 1.
            vsim(1) = vpol(k0_ind,A0_ind); % Firm value in period 1 given k0 and A0.
            ksim(1) = kpol(k0_ind,A0_ind); % Capital choice for period 2 given k0 and A0.
            isim(1) = ipol(k0_ind,A0_ind); % Investment in period 1 given k0 and A0.
            rsim(1) = rpol(k0_ind,A0_ind); % Revenue in period 1 given k0 and A0.
            esim(1) = ipol(k0_ind,A0_ind); % Investment ependiture in period 1 given k0 and A0.
            psim(1) = ppol(k0_ind,A0_ind); % Utility in period 1 given k0 and A0.

            A1_ind = find(rand<=cmat(A0_ind,:)); % Draw productivity for next period.
            A0_ind = A1_ind(1);

            %% Simulate endogenous and exogenous variables.

            for j = 2:T*2 % Time loop.
                kt_ind = find(ksim(j-1)==kgrid); % Capital choice in the previous period is the state today. Find where the latter is on the grid.
                Asim(j) = Agrid(A0_ind); % Productivity in period t.
                vsim(j) = vpol(kt_ind,A0_ind); % Firm value in period t.
                ksim(j) = kpol(kt_ind,A0_ind); % Capital stock for period t+1.
                isim(j) = ipol(kt_ind,A0_ind); % Investment in period t.
                rsim(j) = rpol(kt_ind,A0_ind); % Revenue in period t.
                esim(j) = epol(kt_ind,A0_ind); % Investment expenditure in period t.
                psim(j) = ppol(kt_ind,A0_ind); % Profit in period t.
                A1_ind = find(rand<=cmat(A0_ind,:)); % Draw next state.
                A0_ind = A1_ind(1); % State next period.
            end

            sim = struct();
            
            % Burn the first half.
            sim.Asim = Asim(T+1:2*T,1); % Simulated productivity.
            sim.vsim = vsim(T+1:2*T,1); % Simulated output.
            sim.ksim = ksim(T+1:2*T,1); % Simulated capital choice.
            sim.isim = isim(T+1:2*T,1); % Simulated investment.
            sim.rsim = rsim(T+1:2*T,1); % Simulated revenue.
            sim.esim = esim(T+1:2*T,1); % Simulated investment expenditure.
            sim.psim = psim(T+1:2*T,1); % Simulated profit.
             
        end
        
    end
end