%% File Info.

%{

    simulate.m
    ----------
    This code simulates the model.

%}

%% Simulate class.

classdef q_simulate
    methods(Static)
        %% Simulate the model. 
        
        function sim = lc(par,sol)            
            %% Set up.
            
            agrid = par.agrid; % Assets today (state variable).

            apol = sol.a; % Policy function for capital.
            cpol = sol.c; % Policy function for consumption.

            TT = par.TT; % Time periods.
            NN = par.NN; % People.
            T = par.T; % Life span.
            tr = par.t_r; % Retirement.

            kappa = par.k; % Share of income as pension.
            ygrid = par.ygrid; % Exogenous income.
            pmat = par.pmat; % Transition matrix.

            ysim = nan(TT,NN); % Container for simulated income.
            Asim = nan(TT,NN); % Container for simulated savings.
            tsim = nan(TT,NN); % Container for simulated age.
            csim = nan(TT,NN); % Container for simulated consumption.
            usim = nan(TT,NN); % Container for simulated utility.
            ksim = nan(TT,NN);
            isim = nan(TT,NN);
            
          
            rng(par.seed);

            pmat0 = pmat^100; % Stationary distirbution.
            cmat = cumsum(pmat,2); % CDF matrix.

            y0_ind = randsample(par.ylen,NN,true,pmat0(1,:))'; % Index for initial income.
            a0_ind = randsample(par.alen,NN,true)'; % Index for initial wealth.
            t0_ind = randsample(T,NN,true)'; % Index for initial wealth.
            yr = nan(NN,1); % Retirement income

            for i = 1:NN % Person loop.
                
                if t0_ind(i)>=tr % Retired now.
                    yr(i) = ygrid(y0_ind(i)); % Store for pension.
                    ysim(1,i) = kappa.*yr(i); % Pension in period 0 given age.
                else
                    ysim(1,i) = ygrid(y0_ind(i));
                end

                tsim(1,i) = t0_ind(i); % Age in period 0.
                csim(1,i) = cpol(a0_ind(i),t0_ind(i),y0_ind(i)); % Consumption in period 0 given a0.
                Asim(1,i) = apol(a0_ind(i),t0_ind(i),y0_ind(i)); % Savings for period 1 given a0.

                if t0_ind(i) == tr-1 % Retired next period.
                    yr(i) = ygrid(y0_ind(i)); % Store as pension for next period
                elseif t0_ind(i) < tr-1
                    y1_ind = find(rand<=cmat(y0_ind(i),:)); % Draw income shock for next period.
                    y0_ind(i) = y1_ind(1);
                end

            end

            usim(1,:) = q_model.utility(csim(1,:),par); % Utility in period 0 given a0.

            %% Simulate endogenous variables.

            for j = 2:TT % Time loop.
                for i = 1:NN % Person loop.

                    age = tsim(j-1,i)+1; % Age in period t.

                    if age <= T % Check if still alive.
                        
                        if age>=tr % Retired
                            ysim(j,i) = kappa.*yr(i); % Pension in period t given age.
                        else
                            ysim(j,i) = ygrid(y0_ind(i)); % Pension in period t given age.
                        end

                        tsim(j,i) = age; % Age in period t.
                        at_ind = find(Asim(j-1,i)==agrid); % Savings choice in the previous period is the state today. Find where the latter is on the grid.
                        %kt_ind = find(ksim(j-1,i)==kgrid);
                        csim(j,i) = cpol(at_ind,age,y0_ind(i)); % Consumption in period t.
                        Asim(j,i) = apol(at_ind,age,y0_ind(i)); % Savings for period t+1.
                        %ksim(j) = kpol(kt_ind,A0_ind); % Capital stock for period t+1.
                        %isim(j) = ipol(kt_ind,A0_ind); % Investment in period t.
                        usim(j) = q_model.utility(csim(j),par); % Utility in period t.
                        if age == tr-1 % Retire next period
                            yr(i) = ygrid(y0_ind(i)); % Store as pension for next period
                        elseif age < tr-1
                            y1_ind = find(rand<=cmat(y0_ind(i),:)); % Draw income shock for next period.
                            y0_ind(i) = y1_ind(1);
                        end

                    end
                end
            end

            sim = struct();
            
            sim.ysim = ysim; % Simulated output.
            sim.Asim = Asim; % Simulated savings.
            sim.tsim = tsim; % Simulated age.
            sim.csim = csim; % Simulated consumption.
            sim.usim = usim; % Simulated utility.
            %sim.ksim = ksim;
             
        end
        
    end
end