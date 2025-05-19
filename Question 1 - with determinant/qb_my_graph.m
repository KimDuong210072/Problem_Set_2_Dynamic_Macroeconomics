classdef qb_my_graph
    methods(Static)

        %% Plot policy and value functions
        function [] = plot_policy(par, sol, sim)
            %% Setup
            ystate = par.ygrid;
            age = 1:par.T;

            %% Consumption Policy Function
            figure(1)
            surf(age(1:5:end), ystate, squeeze(sol.c(1,1:5:end,:))')
            xlabel('$t$', 'Interpreter', 'latex')
            ylabel('$y_t$', 'Interpreter', 'latex') 
            zlabel('$c_t$', 'Interpreter', 'latex') 
            title('Consumption Policy Function: Lowest $a_t$', 'Interpreter', 'latex')
            shading interp; view(135,30)

            figure(2)
            surf(age(1:5:end), ystate, squeeze(sol.c(end,1:5:end,:))')
            xlabel('$t$', 'Interpreter', 'latex')
            ylabel('$y_t$', 'Interpreter', 'latex') 
            zlabel('$c_t$', 'Interpreter', 'latex') 
            title('Consumption Policy Function: Highest $a_t$', 'Interpreter', 'latex')
            shading interp; view(135,30)

            %% Saving Policy Function
            figure(3)
            surf(age(1:5:end), ystate, squeeze(sol.a(1,1:5:end,:))')
            xlabel('$t$', 'Interpreter', 'latex')
            ylabel('$y_t$', 'Interpreter', 'latex') 
            zlabel('$a_{t+1}$', 'Interpreter', 'latex') 
            title('Saving Policy Function: Lowest $a_t$', 'Interpreter', 'latex')
            shading interp; view(135,30)

            figure(4)
            surf(age(1:5:end), ystate, squeeze(sol.a(end,1:5:end,:))')
            xlabel('$t$', 'Interpreter', 'latex')
            ylabel('$y_t$', 'Interpreter', 'latex') 
            zlabel('$a_{t+1}$', 'Interpreter', 'latex') 
            title('Saving Policy Function: Highest $a_t$', 'Interpreter', 'latex')
            shading interp; view(135,30)

            %% Value Function
            if isfield(sol, 'v')
                figure(5)
                surf(age(1:5:end), ystate, squeeze(sol.v(1,1:5:end,:))')
                xlabel('$t$', 'Interpreter', 'latex')
                ylabel('$y_t$', 'Interpreter', 'latex') 
                zlabel('$v_t(a_t,t)$', 'Interpreter', 'latex')
                title('Value Function: Lowest $a_t$', 'Interpreter', 'latex')
                shading interp; view(135,30)

                figure(6)
                surf(age(1:5:end), ystate, squeeze(sol.v(end,1:5:end,:))')
                xlabel('$t$', 'Interpreter', 'latex')
                ylabel('$y_t$', 'Interpreter', 'latex') 
                zlabel('$v_t(a_t,t)$', 'Interpreter', 'latex')
                title('Value Function: Highest $a_t$', 'Interpreter', 'latex')
                shading interp; view(135,30)
            end

            %% Life-Cycle Profiles (Smoothed)
            lcp_c = qb_my_graph.mean_by_age(sim.csim, sim.tsim, par.T);
            lcp_a = qb_my_graph.mean_by_age(sim.Asim, sim.tsim, par.T);
            lcp_u = qb_my_graph.mean_by_age(sim.usim, sim.tsim, par.T);
            lcp_g = qb_my_graph.mean_by_age(sim.goldsim, sim.tsim, par.T);

            figure(7)
            plot(age, lcp_c, 'LineWidth', 1.5)
            xlabel('Age', 'Interpreter', 'latex')
            ylabel('$c^{\mathrm{sim}}_t$', 'Interpreter', 'latex')
            title('LCP of Consumption')
            grid on
            
            figure(8)
            plot(age, lcp_a, 'b', 'LineWidth', 1.5)
            xlabel('Age', 'Interpreter', 'latex')
            ylabel('$a^{\mathrm{sim}}_{t+1}$', 'Interpreter', 'latex')
            title('LCP of Savings')
            grid on
            
            figure(9)
            plot(age, lcp_u, 'LineWidth', 1.5)
            xlabel('Age', 'Interpreter', 'latex')
            ylabel('$u^{\mathrm{sim}}_t$', 'Interpreter', 'latex')
            title('LCP of Utility')
            grid on
            
            figure(10)
            plot(age, lcp_g, 'LineWidth', 1.5)
            xlabel('Age', 'Interpreter', 'latex')
            ylabel('$g^{\mathrm{sim}}_t$', 'Interpreter', 'latex')
            title('LCP of Gold Holdings')
            grid on
        end

        %% Smoothed mean by age with minimum sample count
        function profile = mean_by_age(data, tsim, T)
            profile = nan(T, 1);
            valid = ~isnan(tsim) & tsim >= 1 & tsim <= T;
            min_obs = 30;  % Require at least 30 people at each age

            for t = 1:T
                mask = valid & (tsim == t) & ~isnan(data);
                if sum(mask(:)) >= min_obs
                    profile(t) = mean(data(mask), 'omitnan');
                end
            end

            % Smooth the result using moving average
            profile = movmean(profile, 3, 'omitnan');
        end
    end
end
