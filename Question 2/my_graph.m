%% File Info.
%{
    my_graph.m
    ----------
    This code plots the value and policy functions and the time path of the variables.
%}

%% Graph class.

classdef my_graph
    methods(Static)
        function [] = plot_policy(par, sol, sim)
            %% Fix p index to plot mid-price slice
            p_index = ceil(par.plen / 2);
            p_val = par.pgrid(p_index);
            
            %% Plot capital policy function across A
            figure(1)
            hold on
            for A_index = 1:par.Alen
                plot(par.kgrid, squeeze(sol.k(:, A_index, p_index)), ...
                     'DisplayName', ['$A = $', num2str(par.Agrid(A_index), '%.2f')])
            end
            xlabel({'$k_t$'}, 'Interpreter', 'latex')
            ylabel({'$k_{t+1}$'}, 'Interpreter', 'latex')
            title('Capital Policy Function across Productivity States', 'Interpreter', 'latex')
            legend('Interpreter','latex', 'Location','best')
            hold off
            
            %% Plot investment policy function across A
            figure(2)
            hold on
            for A_index = 1:par.Alen
                plot(par.kgrid, squeeze(sol.i(:, A_index, p_index)), ...
                     'DisplayName', ['$A = $', num2str(par.Agrid(A_index), '%.2f')])
            end
            xlabel({'$k_t$'}, 'Interpreter', 'latex')
            ylabel({'$i_{t}$'}, 'Interpreter', 'latex')
            title('Investment Policy Function across Productivity States', 'Interpreter', 'latex')
            legend('Interpreter','latex', 'Location','best')
            hold off
            
            %% Plot revenue function across A
            figure(3)
            hold on
            for A_index = 1:par.Alen
                plot(par.kgrid, squeeze(sol.r(:, A_index, p_index)), ...
                     'DisplayName', ['$A = $', num2str(par.Agrid(A_index), '%.2f')])
            end
            xlabel({'$k_t$'}, 'Interpreter', 'latex')
            ylabel({'$r_t$'}, 'Interpreter', 'latex')
            title('Revenue Function across Productivity States', 'Interpreter', 'latex')
            legend('Interpreter','latex', 'Location','best')
            hold off
            
            %% Plot expenditure function across A
            figure(4)
            hold on
            for A_index = 1:par.Alen
                plot(par.kgrid, squeeze(sol.e(:, A_index, p_index)), ...
                     'DisplayName', ['$A = $', num2str(par.Agrid(A_index), '%.2f')])
            end
            xlabel({'$k_t$'}, 'Interpreter', 'latex')
            ylabel({'$e_t$'}, 'Interpreter', 'latex')
            title('Expenditure Function across Productivity States', 'Interpreter', 'latex')
            legend('Interpreter', 'latex', 'Location', 'best')
            hold off
            
            %% Plot profit function across A
            figure(5)
            hold on
            for A_index = 1:par.Alen
                plot(par.kgrid, squeeze(sol.p(:, A_index, p_index)), ...
                     'DisplayName', ['$A = $', num2str(par.Agrid(A_index), '%.2f')])
            end
            xlabel({'$k_t$'}, 'Interpreter', 'latex')
            ylabel({'$\pi_t$'}, 'Interpreter', 'latex')
            title('Profit Function across Productivity States', 'Interpreter', 'latex')
            legend('Interpreter','latex', 'Location','best')
            hold off
            
            %% Plot value function across A
            figure(6)
            hold on
            for A_index = 1:par.Alen
                plot(par.kgrid, squeeze(sol.v(:, A_index, p_index)), ...
                     'DisplayName', ['$A = $', num2str(par.Agrid(A_index), '%.2f')])
            end
            xlabel({'$k_t$'}, 'Interpreter', 'latex')
            ylabel({'$v_t$'}, 'Interpreter', 'latex')
            title('Value Function across Productivity States', 'Interpreter', 'latex')
            legend('Interpreter','latex', 'Location','best')
            hold off

            %% Time vector after burn-in
            T = par.T;
            tgrid = 1:T;

            %% Time path plots (averaged over firms)

            % Simulated productivity
            figure(7)
            plot(tgrid, mean(sim.Asim, 2))
            xlabel({'Time'}, 'Interpreter', 'latex')
            ylabel({'$A_t^{sim}$'}, 'Interpreter', 'latex')
            title('Average Simulated Productivity States')

            % Simulated prices
            figure(8)
            plot(tgrid, mean(sim.Psim, 2))
            xlabel({'Time'}, 'Interpreter', 'latex')
            ylabel({'$p_t^{sim}$'}, 'Interpreter', 'latex')
            title('Average Simulated Price Shocks')

            % Simulated capital
            figure(9)
            plot(tgrid, mean(sim.ksim, 2))
            xlabel({'Time'}, 'Interpreter', 'latex')
            ylabel({'$k_t^{sim}$'}, 'Interpreter', 'latex')
            title('Average Simulated Capital')

            % Simulated investment
            figure(10)
            plot(tgrid, mean(sim.isim, 2))
            xlabel({'Time'}, 'Interpreter', 'latex')
            ylabel({'$i_t^{sim}$'}, 'Interpreter', 'latex')
            title('Average Simulated Investment')

            % Simulated expenditure
            figure(11)
            plot(tgrid, mean(sim.esim, 2))
            xlabel({'Time'}, 'Interpreter', 'latex')
            ylabel({'$e_t^{sim}$'}, 'Interpreter', 'latex')
            title('Average Simulated Expenditure')

            % Simulated revenue
            figure(12)
            plot(tgrid, mean(sim.rsim, 2))
            xlabel({'Time'}, 'Interpreter', 'latex')
            ylabel({'$r_t^{sim}$'}, 'Interpreter', 'latex')
            title('Average Simulated Revenue')

            % Simulated profit
            figure(13)
            plot(tgrid, mean(sim.psim, 2))
            xlabel({'Time'}, 'Interpreter', 'latex')
            ylabel({'$\pi_t^{sim}$'}, 'Interpreter', 'latex')
            title('Average Simulated Profit')

            % Simulated firm value
            figure(14)
            plot(tgrid, mean(sim.vsim, 2))
            xlabel({'Time'}, 'Interpreter', 'latex')
            ylabel({'$v_t^{sim}$'}, 'Interpreter', 'latex')
            title('Average Simulated Firm Value')
        end
    end
end
