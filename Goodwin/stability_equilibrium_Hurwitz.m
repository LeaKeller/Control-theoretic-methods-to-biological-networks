function SE = stability_equilibriun_Hurwitz()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the function of det(H_2) for different values of
% n. H_2 represents the Hurwitz matrix. For the Hurwitz criterion to be
% satisfied and obtain stability, det(H_2)>0.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    syms z

    k1 = 1;
    k2 = 0.1;
    k3 = 1;
    k4 = 0.1;
    k5 = 1;
    k6 = 0.1;
    ki = 1;

    a = -k2/(k1*k3*k5/ki)^(1/3); 
    b = -k4/(k1*k3*k5/ki)^(1/3);
    c = -k6/(k1*k3*k5/ki)^(1/3);
    
    color = ['r','g','b','c','n','y','k','w',[0.6350 0.0780 0.1840],[0.3010 0.7450 0.9330]];
    
    figure;
    for n = 1:10
        f = @(z) (a+b+c)*(a*b+b*c+c*a) - c*a*b - (n*z^(n-1))/((1+z^n)^2); hold on;
        fplot(f);
    end
    yline(0,'--r');
    legend('n=1','n=2','n=3','n=4','n=5','n=6','n=7','n=8','n=9','n=10','q_1q_2-q_3=0');
    xlim([-2 0]);
    ylim([-20 20]);
end