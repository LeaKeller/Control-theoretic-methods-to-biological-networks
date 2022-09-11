function H = Hopf_Brusselator()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the Hopf bifurcation of the Brusselator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % compute the Jacobian matrix at the equilibrium point
    syms a b
    Jac = [b-1 a^2 ; -b -a^2];

    % compute the eigenvalues of the Jacobian 
    % First compute the characteristic polynomial 
    syms lam
    Jac_l = [b-1-lam a^2 ; -b -a^2-lam];
    lamP = charpoly(Jac,lam);
    eigenvalues = eig(Jac_l);

    format rat

    n = 14;
    a=3; 
    b = [];
    for j = 1:n
        b = [b, 1.5*j];
    end
    c = ['r','g','b','c','m','y','k','w',[0.6350 0.0780 0.1840],[0.3010 0.7450 0.9330],[0.4660 0.6740 0.1880],[0.4940 0.1840 0.5560],[0.9290 0.6940 0.1250],[0.8500 0.3250 0.0980],[0 0.4470 0.7410]];
    
    figure;
    for i = 1:n
        txt = ['b =', num2str(b(i))];
        l1(i) = b(i)/2 - (-(- a^2 + 2*a + b(i) - 1)*(a^2 + 2*a - b(i) + 1))^(1/2)/2 - a^2/2 - 1/2;
        l2(i) = b(i)/2 + (-(- a^2 + 2*a + b(i) - 1)*(a^2 + 2*a - b(i) + 1))^(1/2)/2 - a^2/2 - 1/2;
        scatter([real(l1(i)),real(l2(i))], [imag(l1(i)),imag(l2(i))], c(i), 'DisplayName', txt); hold on;
    end
    title('a=3');
    plot(real(l1),imag(l1),'DisplayName','\lambda_1');
    plot(real(l2),imag(l2),'DisplayName','\lambda_2');
    xlabel('\Re(\lambda)');
    ylabel('\Im(\lambda)');
    legend show;
end