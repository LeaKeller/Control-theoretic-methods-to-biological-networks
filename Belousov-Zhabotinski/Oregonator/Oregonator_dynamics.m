function OD = Oregonator_dynamics()

    %% Values taken from http://www.scholarpedia.org/article/Oregonator
    %% Case of the system from the website
    
    %% Method
    % ode45 (first try) : too much noise since ode45 is
    % a solver used to solve non stiff differential equations
    % other solvers : example ode23s (best) : solve stiff differntial
    % equations
    %%
    
    A = 0.06;
    B = 0.02;
    
    k1 = 1.28;
    k2 = 2.4*10^6;
    k3 = 33.6;
    k4 = 2400;
    k5 = 1;
    
    f = [0.1, 1, 1+sqrt(2)];
    
    tfinal = 1500;
    
    % Inintial conditions for full (non-reduced system)
    x01 = 1;
    x02 = 1;
    x03 = 1;
    
    for i=1:length(f)
        [t, x] = ode23s(@(t, x) BZ_Oregonator_system(t,x,k1,k2,k3,k4,k5,A,B,f(i)) , [0 tfinal], [x01 ; x02 ; x03]);

        x1 = x(:,1);
        x2 = x(:,2);
        x3 = x(:,3);

        tau = k5*B*t; % scaled t

        figure;
        plot(tau,log10(x1));
        xlabel('\tau');
        ylabel('log_{10}(x)');
        title(['f = ', num2str(f(i))]);
        grid on; 

        figure;
        plot(tau,log10(x2));
        xlabel('\tau');
        ylabel('log_{10}(y)');
        title(['f = ', num2str(f(i))]);
        grid on; 

        figure;
        plot(tau,log10(x3));
        xlabel('\tau');
        ylabel('log_{10}(z)');
        title(['f = ', num2str(f(i))]);
        grid on;
    end
end