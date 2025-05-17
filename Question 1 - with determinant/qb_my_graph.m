%% File Info.

%{

    my_graph.m
    ----------
    This code plots the value and policy functions and the time path of the variables.

%}

%% Graph class.

classdef qb_my_graph
    methods(Static)
        %% Plot value and policy functions.
        
        function [] = plot_policy(par,sol,sim)
            %% Plot consumption policy function.

            ystate = par.ygrid;
            age = linspace(1,par.T,par.T);
            
            figure(1)
            
            surf(age(1:5:end),ystate,squeeze(sol.c(1,1:5:end,:))')
                xlabel({'t'},'Interpreter','latex')
                ylabel({'$y_{t}$'},'Interpreter','latex') 
                zlabel({'$c_{t}$'},'Interpreter','latex') 
            title('Consumption Policy Function, Lowest $a_t$','Interpreter','latex')
            
            figure(2)
            
            surf(age(1:5:end),ystate,squeeze(sol.c(end,1:5:end,:))')
                xlabel({'t'},'Interpreter','latex')
                ylabel({'$y_{t}$'},'Interpreter','latex') 
                zlabel({'$c_{t}$'},'Interpreter','latex') 
            title('Consumption Policy Function, Highest $a_t$','Interpreter','latex')
            
            %% Plot saving policy function.
            
            figure(3)
            
            surf(age(1:5:end),ystate,squeeze(sol.a(1,1:5:end,:))')
                xlabel({'t'},'Interpreter','latex')
                ylabel({'$y_{t}$'},'Interpreter','latex') 
                zlabel({'$a_{t+1}$'},'Interpreter','latex') 
            title('Saving Policy Function, Lowest $a_t$','Interpreter','latex')
            
            figure(4)
            
            surf(age(1:5:end),ystate,squeeze(sol.a(end,1:5:end,:))')
                xlabel({'t'},'Interpreter','latex')
                ylabel({'$y_{t}$'},'Interpreter','latex') 
                zlabel({'$a_{t+1}$'},'Interpreter','latex') 
            title('Saving Policy Function, Highest $a_t$','Interpreter','latex')
            
            %% Plot value function.
            
            figure(5)
            
            surf(age(1:5:end),ystate,squeeze(sol.v(1,1:5:end,:))')
                xlabel({'t'},'Interpreter','latex')
                ylabel({'$y_{t}$'},'Interpreter','latex') 
                zlabel({'$v_t(a_t,t)$'},'Interpreter','latex')
            title('Value Function, Lowest $a_t$','Interpreter','latex')

            figure(6)
            
            surf(age(1:5:end),ystate,squeeze(sol.v(end,1:5:end,:))')
                xlabel({'t'},'Interpreter','latex')
                ylabel({'$y_{t}$'},'Interpreter','latex') 
                zlabel({'$v_t(a_t,t)$'},'Interpreter','latex')
            title('Value Function, Highest $a_t$','Interpreter','latex')

            %% Plot consumption policy function.
            
            lcp_c = nan(par.T,1);
            lcp_a = nan(par.T,1);
            lcp_u = nan(par.T,1);
            
                for i=1:par.T
                    lcp_c(i) = mean(sim.csim(sim.tsim==i),"omitnan");
                    lcp_a(i) = mean(sim.Asim(sim.tsim==i),"omitnan");
                    lcp_u(i) = mean(sim.usim(sim.tsim==i),"omitnan");
                end
            figure(7)
            plot(age,lcp_c)
                    xlabel({'$Age$'},'Interpreter','latex')
                    ylabel({'$c^{sim}_{t}$'},'Interpreter','latex')
            title('LCP of Consumption')
            %% Plot saving policy function.
            figure(8)
            plot(age,lcp_a)
                    xlabel({'$Age$'},'Interpreter','latex')
                    ylabel({'$a^{sim}_{t+1}$'},'Interpreter','latex')
            title('LCP of Savings')
            %% Plot value function.
            figure(9)
            plot(age,lcp_u)
                    xlabel({'$Age$'},'Interpreter','latex')
                    ylabel({'$u^{sim}_t$'},'Interpreter','latex')
            title('LCP of Utility')
            % Plot consumption profiles
                %for t = 1:par.T
                    %consumption_profiles(t) = mean(sim.csim(sim.tsim == t), 'omitnan');
                    %wealth_profiles(t) = mean(sim.Asim(sim.tsim == t), 'omitnan');
                %end
                %figure;
                %hold on;
                %for i = 1:length(beta_values)
                    %plot(1:T, consumption_profiles(i,:), 'DisplayName', ['\beta = ' num2str(beta_values(i))]);
                %end
                %xlabel('Age');
                %ylabel('Average Consumption');
                %title('Life-Cycle Consumption Profiles for \gamma = 2.00');
                %legend show;
                %grid on;
                %hold off;
                
                % Plot wealth profiles
                %figure;
                %hold on;
                %for i = 1:length(beta_values)
                    %plot(1:T, wealth_profiles(i,:), 'DisplayName', ['\beta = ' num2str(beta_values(i))]);
                %end
                %xlabel('Age');
                %ylabel('Average Wealth');
                %title('Life-Cycle Wealth Profiles for \gamma = 2.00');
                %legend show;
                %grid on;
                %hold off;
            end
      end
end
