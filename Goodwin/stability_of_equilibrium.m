function S = stability_of_equilibrium()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function computing the stability of the equilibria. This is done by
% checking if the real parts of the eigenvalues are in the unit circle. An
% other plot shows the dynamics of the system to the equilibrium.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    syms x y z

    % Parameters
    k1 = 1;
    k2 = 0.1;
    k3 = 1;
    k4 = 0.1;
    k5 = 1;
    k6 = 0.1;
    ki = 1;

    n = 10;

    tfinal = 500;
    
    x0 = 1;
    y0 = 1;
    z0 = 1;
    
    [t, X] = ode23s(@(t, X) Goodwin(t,X,k1,k2,k3,k4,k5,k6,ki,n), [0 tfinal], [x0 ; y0 ; z0]);
    
    x1 = X(:,1);
    x2 = X(:,2);
    x3 = X(:,3);
    
    e1 = k1*(ki^n/(ki^n+z^n)) - k2*x;
    e2 = k3*x - k4*y;
    e3 = k5*y - k6*z;
    
    solution = vpasolve([e1 == 0, e2 == 0, e3 == 0], [x y z]); % equilibrium points
    
    latextable = array2table([solution.x, solution.y, solution.z], 'VariableNames',{'X','Y','Z'});
    
    figure;
    plot3(x1,x2,x3); hold on;
    for i=1:length(size(solution.x,1))
        scatter3(solution.x(i),solution.y(i),solution.z(i)); hold on;
    end
    grid on;
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    legend('Dynamics','Equilibria');
    
    figure;
    for i=1:size(solution.x,1)
        scatter([real(solution.x(i)),real(solution.y(i)),real(solution.z(i))], [imag(solution.x(i)),imag(solution.y(i)),imag(solution.z(i))]); hold on;
    end
    
    list_eig1_real = [];
    list_eig2_real = [];
    list_eig3_real = [];
    list_eig1_im = [];
    list_eig2_im = [];
    list_eig3_im = [];
    
    for i = 1:size(solution.x,1)
        n = i;
        J = [-k2 0 ki*(ki^n*n*solution.z(i)^(n-1))/(ki^n+solution.z(i)^n)^2 ; k3 -k4 0 ; 0 k5 -k6]; % Jacobian matrix for the first solution couples
        eigenvalues = eig(J);
        list_eig1_real = [list_eig1_real real(eigenvalues(1))];
        list_eig2_real = [list_eig2_real real(eigenvalues(2))];
        list_eig3_real = [list_eig3_real real(eigenvalues(3))];
        list_eig1_im = [list_eig1_im imag(eigenvalues(1))];
        list_eig2_im = [list_eig2_im imag(eigenvalues(2))];
        list_eig3_im = [list_eig3_im imag(eigenvalues(3))];
        fprintf('First real part of the eigenvalues corresponding to the 11 different equilibrium points = %g\n', eigenvalues(1));
        fprintf('Second real part of the eigenvalues corresponding to the 11 different equilibrium points = %g\n', eigenvalues(2));
        fprintf('Third real part of the eigenvalues corresponding to the 11 different equilibrium points = %g\n', eigenvalues(3));
    end
    
    L1_real = list_eig1_real;
    L2_real = list_eig2_real;
    L3_real = list_eig3_real;
    L1_im = list_eig1_im;
    L2_im = list_eig2_im;
    L3_im = list_eig3_im;

    figure;
    for i = 1:size(list_eig1_real,1)*size(list_eig1_real,2)
        scatter([L1_real(i), L2_real(i), L3_real(i)], [L1_im(i), L2_im(i), L3_im(i)]); hold on;
    end
    circle = @(X,Y) X.^2 + Y.^2 -1;
    ezplot(circle); hold on; % unit circle
    axis([-1 1 -1 1]);
       
   for i = 1:size(solution.x,1)
       n = i;
       J = [-k2 0 ki*(ki^n*n*solution.z(i)^(n-1))/(ki^n+solution.z(i)^n)^2 ; k3 -k4 0 ; 0 k5 -k6]; % Jacobian matrix for the first solution couples
       tr = trace(J);
       fprintf('tr(J) = %g\n', tr);
    end

end