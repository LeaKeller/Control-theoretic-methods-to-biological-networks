function S = Simulate_controled_system

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function plots the dynamics of the Brusselator by adding a
% controlers to the concentrations x and y. This is performed using LQR.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    a = 3;
    b = 16.5; % choose b such that the equilibrium is unstable

    syms x y ux uy k1 k2
    xdot = a - (b+1)*(x+ux) + (x+ux)^2*(y+uy);
    ydot = b*(x+ux) - (x+ux)^2*(y+uy);

    Xdot = [xdot;ydot];
    J = jacobian(Xdot,[x y ux uy]);

    J_equilibrium = subs(J,[x y ux uy],[a b/a 0 0]);

    A = J_equilibrium(1:2,1:2);
    B = J_equilibrium(1:2,3:4);

    K = [k1 0;0 k1];

    M = A-B*K;
    eigenval_M = eig(M);
    [V,D,W] = eig(M);

    figure;
    fplot(eigenval_M);
    hold on;
    yline(0,'--','y=0');
    xline(1,'--','x=1');
    legend('\lambda_1','\lambda_2');
    grid on;

    A_double = double(A);
    B_double = double(B);
    state_model = ss(A_double, B_double, eye(2), zeros(2));
    Q = eye(2)*1000;
    R = eye(2)*1;
    K_lqr = lqr(state_model,Q,R);

    % Figures
    
    tfinal = 50;
    K = eye(2)*1;

    x_eq = [a; b/a];
    %x_0 = x_eq;
    x_0 = x_eq - 0.1;
    options = odeset('AbsTol',1e-10,'RelTol',1e-10);

    [t, x] = ode45(@(t, x) BZ_reduced_test(t, x, a, b, K_lqr) , [0 tfinal], x_0, options);%[x01 ; x02]);

    x1 = x(:,1);
    x2 = x(:,2);

    figure;
    plot(x1,x2); hold on;
    scatter(x_0(1),x_0(2),'filled');
    scatter(x_eq(1),x_eq(2),'filled');
    legend('Reaction mouvement','Starting point','Equilibrium');
    xlabel('x_1');
    ylabel('x_2');
    grid on;
    title(['a=3, ','b=', num2str(b)]);

    figure;
    plot(t,x1); hold on; 
    plot(t,x2);
    xlabel('t');
    legend('x_1','x_2');
    grid on; 
    title(['a=3, ','b=', num2str(b)]);a = 3;
    b = 16.5;

    syms x y ux uy k1 k2
    xdot = a - (b+1)*(x+ux) + (x+ux)^2*(y+uy);
    ydot = b*(x+ux) - (x+ux)^2*(y+uy);

    Xdot = [xdot;ydot];
    J = jacobian(Xdot,[x y ux uy]);

    J_equilibrium = subs(J,[x y ux uy],[a b/a 0 0]);

    A = J_equilibrium(1:2,1:2);
    B = J_equilibrium(1:2,3:4);

    K = [k1 0;0 k1];

    M = A-B*K;
    eigenval_M = eig(M);
    [V,D,W] = eig(M);
end