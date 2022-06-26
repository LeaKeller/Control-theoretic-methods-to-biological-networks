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
    
    f = 1;
    
    tfinal = 1500;
    
    x01 = 1;
    x02 = 1;
    x03 = 1;
    
    [t, x] = ode23s(@(t, x) BZ_Oregonator_reduced(t,x,k1,k2,k3,k4,k5,A,B,f) , [0 tfinal], [x01 ; x02 ; x03]);
    
%     x1 = (k3*A*x(:,1))/(2*k4); % scaled x
%     x2 = (k3*A*x(:,2))/k2; % scaled y
%     x3 = (k3*A)^2*x(:,3)/(k5*k4*B); % scaled z
    
    x1 = x(:,1);
    x2 = x(:,2);
    x3 = x(:,3);
    
    tau = k5*B*t; % scaled t
    
    figure;
    plot(tau,log10(x1));
    xlabel('\tau');
    ylabel('log_{10}(x)');
    grid on; 
    
    figure;
    plot(tau,log10(x2));
    xlabel('\tau');
    ylabel('log_{10}(y)');
    grid on; 
    
    figure;
    plot(tau,log10(x3));
    xlabel('\tau');
    ylabel('log_{10}(z)');
    grid on;
end