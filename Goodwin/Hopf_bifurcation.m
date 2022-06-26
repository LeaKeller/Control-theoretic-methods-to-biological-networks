%% Study Hopf bifurcation

function H = Hopf_bifurcation()
    syms x y z 
    
    % Parameters
    k1 = 1;
    k2 = 0.1;
    k3 = 1;
    k4 = 0.1;
    k5 = 1;
    k6 = 0.1;
    ki = 1;
    
    figure;
    for j = 1:14
    n = j;

    e1 = k1*(ki^n/(ki^n+z^n)) - k2*x;
    e2 = k3*x - k4*y;
    e3 = k5*y - k6*z;

    solution = vpasolve([e1 == 0, e2 == 0, e3 == 0], [x y z]);

    c = ['r','g','b','c','m','y','k',[0.6350 0.0780 0.1840],[0.3010 0.7450 0.9330],[0.4660 0.6740 0.1880],[0.4940 0.1840 0.5560],[0.9290 0.6940 0.1250],[0.8500 0.3250 0.0980],[0 0.4470 0.7410]];

    % Hopf bifurcation diagram for the 3 different eigenvalues for different values
    % of n
        for i = 1:size(solution.x,1)
            txt = ['n =', num2str(j)];
            J = [-k2 0 ki*(ki^n*n*solution.z(i)^(n-1))/(ki^n+solution.z(i)^n)^2 ; k3 -k4 0 ; 0 k5 -k6]; % Jacobian matrix
            eigenvalues = eig(J);
         end
        scatter([real(eigenvalues(1)), real(eigenvalues(2)), real(eigenvalues(3))], [imag(eigenvalues(1)), imag(eigenvalues(2)), imag(eigenvalues(3))], c(j), 'DisplayName', txt); hold on;
    end
    xlabel('Real(\lambda)');
    ylabel('Im(\lambda)');
    grid on;
    legend show;
end