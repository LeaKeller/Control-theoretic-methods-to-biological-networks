function GD = Gyorgyri_dynamics()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Concentrations and rate constants values taken from 
% REVISED DYNAMICS OF THE BELOUSOV-ZHABOTINSKY REACTION MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    A = 0.1;
    H = 0.26;
    M = 0.25;
    C = 0.000833;
    
    a = 666.7;
    b = 0.3478;
    
    k1 = 4.0e6*M^(-2);
    k2 = 2.0*M^(-3);
    k3 = 3000*M^(-1);
    k4 = 55.2*M^(-5/2);
    k5 = 7000*M^(-1);
    k6 = 0.09*M^(-1);
    k7 = 0.23*M^(-1);
    
    tfinal = 100;
    
    x01 = 1;
    x02 = 1;
    x03 = 1;
    x04 = 1;
    
    %opts = odeset('NormControl','on');%odeset('RelTol',2.371986e-20,'AbsTol',2.371986e-20);
    [t, x] = ode45(@(t, x) BZ_Gyorgyri_system(t,x,k1,k2,k3,k4,k5,k6,k7,A,H,C,M,a,b), [0 tfinal], [x01 ; x02 ; x03 ; x04]);
    
    x1 = x(:,1);
    x2 = x(:,2);
    x3 = x(:,3);
    x4 = x(:,4);
    
    figure;
    plot(t,log(x1));
    xlabel('t');
    ylabel('log(x)');
    grid on; 
    
    figure;
    plot(t,log(x2));
    xlabel('t');
    ylabel('log(y)');
    grid on; 
    
    figure;
    plot(t,log(x3)); 
    xlabel('t');
    ylabel('log(z)');
    grid on;
    
    figure;
    plot(t,log(x4)); 
    xlabel('t');
    ylabel('log(v)');
    grid on;
    
    % Trajectories  
    figure;
    plot3(x1,x2,x4);
    xlabel('x');
    ylabel('z');
    zlabel('v');
    grid on;
end