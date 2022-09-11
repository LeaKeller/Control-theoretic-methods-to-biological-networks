function Stab = Stability_diagram_Hopf()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the stability diagram for the Hopf bifurcation by
% ploting the limit at which this bifurcation occurs with n as a variable.
% This enables to show for which value of n the stability occurs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    syms n
    a = (8/(n-8))^(1/n)*(n/(n-8));
    figure;
    fplot(a, [7 10]);
    ylim([0 100]);
    xlabel('n');
    ylabel('\alpha');
    title('Stability diagram');
    txt_stable = 'Stable';
    txt_unstable = 'Unstable';
    text(7.2,60,txt_stable);
    text(8.9,60,txt_unstable);
    legend('Hopf bifurcation');
end