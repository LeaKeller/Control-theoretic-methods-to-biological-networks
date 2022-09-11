function xdot = BZ_Oregonator_system(t,x,k1,k2,k3,k4,k5,A,B,f)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This function computes the Oregonator system. It is used to plot the 
    % dynamics of the Oregonator in the "Oregonator_dynamics" function.
    % Input : time (t), concentation (x), rate constants (k1, k2, k3, k4, k5),
    % constant concentrations (A and B) and constant (f)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    xdot = [k1*A*x(2) - k2*x(1)*x(2) + k3*A*x(1) - 2*k4*x(1)^2; -k1*A*x(2) - k2*x(1)*x(2) + 1/2*k5*f*B*x(3); 2*k3*A*x(1) - k5*B*x(3)];

%     % 2 variables reduced system
%     q = 2*k1*k4/(k2*k3);
%     eps = k5/(k3*A);
%     xdot = [q/eps*(f*x(2)/(q+x(1))) - 1/eps*(f*x(1)^2/(q+x(1))) + x(1)*(1-x(1)) ; x(1)-x(2)];
end