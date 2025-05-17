%% File Info.

%{

    my_graph.m
    ----------
    This code plots the value and policy functions and the time path of the variables.

%}

%% Graph class.

classdef my_graph
    methods(Static)
        %% Plot value and policy functions.
        
        function [] = plot_policy(par,sol,sim)
            %% Plot capital policy function.

            figure(1)
            
            plot(par.kgrid,sol.k)
                xlabel({'$k_t$'},'Interpreter','latex')
                ylabel({'$k_{t+1}$'},'Interpreter','latex') 
            title('Capital Policy Function','Interpreter','latex')
            
            %% Plot investment policy function.
            
            figure(2)
            
            plot(par.kgrid,sol.i)
                xlabel({'$k_t$'},'Interpreter','latex')
                ylabel({'$i_{t}$'},'Interpreter','latex') 
            title('Investment Policy Function','Interpreter','latex')
            
            %% Plot revenue function.
            
            figure(3)
            
            plot(par.kgrid,sol.r)
                xlabel({'$k_t$'},'Interpreter','latex')
                ylabel({'$r_{t}$'},'Interpreter','latex') 
            title('Revenue Function','Interpreter','latex')
            
            %% Plot expenditure function.
            
            figure(4)
            
            plot(par.kgrid,sol.e)
                xlabel({'$k_t$'},'Interpreter','latex')
                ylabel({'$C(k_{t+1},A_t,k_t)+pI_t$'},'Interpreter','latex') 
            title('Expenditure Function','Interpreter','latex')
            
            %% Plot profit function.
            
            figure(5)
            
            plot(par.kgrid,sol.p)
                xlabel({'$k_t$'},'Interpreter','latex')
                ylabel({'$C(k_{t+1},A_t,k_t)+pi_t$'},'Interpreter','latex') 
            title('Profit Function','Interpreter','latex')
            
            %% Plot value function.
            
            figure(6)
            
            plot(par.kgrid,sol.v)
                xlabel({'$k_t$'},'Interpreter','latex')
                ylabel({'$v_{t}$'},'Interpreter','latex') 
            title('Value Function','Interpreter','latex')
            
            %% Plot simulated revenue shocks.

            tgrid = linspace(1,par.T,par.T);

            figure(7)

            plot(tgrid,sim.Asim)
                xlabel({'Time'},'Interpreter','latex')
                ylabel({'$A^{sim}_t$'},'Interpreter','latex') 
            title('Simulated Revenue Shocks')

            %% Plot simulated capital choice.

            figure(8)

            plot(tgrid,sim.ksim)
                xlabel({'Time'},'Interpreter','latex')
                ylabel({'$k^{sim}_t$'},'Interpreter','latex') 
            title('Simulated Capital Choice')

            %% Plot simulated investment expenditure.

            figure(9)

            plot(tgrid,sim.esim)
                xlabel({'Time'},'Interpreter','latex')
                ylabel({'$C(k^{sim}_{t+1},A^{sim}_t,k^{sim}_t)+pi^{sim}_t$'},'Interpreter','latex') 
            title('Simulated Investment Expenditure')

            %% Plot simulated investment.

            figure(9)

            plot(tgrid,sim.isim)
                xlabel({'Time'},'Interpreter','latex')
                ylabel({'$i^{sim}_t$'},'Interpreter','latex') 
            title('Simulated Investment')

            %% Plot simulated revenue.

            figure(10)

            plot(tgrid,sim.rsim)
                xlabel({'Time'},'Interpreter','latex')
                ylabel({'$y^{sim}_t$'},'Interpreter','latex') 
            title('Simulated Revenue')

            %% Plot simulated profit.

            figure(11)

            plot(tgrid,sim.psim)
                xlabel({'Time'},'Interpreter','latex')
                ylabel({'$\pi^{sim}_t$'},'Interpreter','latex') 
            title('Simulated Profit')

            %% Plot simulated value function.

            figure(12)

            plot(tgrid,sim.vsim)
                xlabel({'Time'},'Interpreter','latex')
                ylabel({'$v^{sim}_t$'},'Interpreter','latex') 
            title('Simulated Firm Value')

        end
        
    end
end