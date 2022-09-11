function D = Dynamics_plot()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function plots the dynamics of the Goodwin model using ode45. Those
% plots are reproduced from the paper "The Goodwin model : behind the Hill
% function".
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Parameters
    k1 = 1;
    k2 = 0.1;
    k3 = 1;
    k4 = 0.1;
    k5 = 1;
    k6 = 0.1;
    ki = 1;

    n = 100; %10;

    % Final time
    tfinal = 200;

    % Initial conditions
    y01 = 0;    
    y02 = 0.2;
    y03 = 2.5;

    [t,y] = ode45(@(t, y) system(t, y, k1, k2, k3, k4, k5, k6, ki, n) , [0 tfinal], [y01; y02; y03]); % Solving the system

    y1 = y(:,1);    % y1 and y2 and y3
    y2 = y(:,2);
    y3 = y(:,3);

    figure;
    plot(t,y1)      % y1 as a function of time (y(t) for the original ODE)
    xlabel('time');
    ylabel('X(t)');
    grid on;

    figure;
    plot(t,y2)      % y2 as a function of time (dy(t)/dt forthe original ODE)
    xlabel('time');
    ylabel('Y(t)');
    grid on;

    figure;
    plot(t,y1)      % y3 as a function of time (z(t) for the original ODE)
    xlabel('time');
    ylabel('Z(t)');
    grid on;

    figure;
    plot(y1,y2)     % Phase plot for y1 and y2
    xlabel('X(t)');
    ylabel('Y(t)');
    grid on;

    figure;
    plot(y2,y3)     % Phase plot for y2 and y3
    xlabel('Y(t)');
    ylabel('Z(t)');
    grid on;

    figure;
    plot(y1,y3)     % Phase plot for y1 and y3
    xlabel('X(t)');
    ylabel('Z(t)');
    grid on;

    figure;
    plot(t,y1); hold on;
    plot(t,y2); hold on;
    plot(t,y3); hold on;
    yline(1,'--r',{'Hill threshold K_i=1'});
    grid on;
    xlabel('time');
    legend('X','Y','Z');
end