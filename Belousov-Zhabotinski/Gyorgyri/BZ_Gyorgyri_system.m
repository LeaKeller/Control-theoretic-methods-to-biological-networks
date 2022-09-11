function xdot = BZ_Gyorgyri_system(t,x,k1,k2,k3,k4,k5,k6,k7,A,H,C,M,a,b)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This function computes the Gyorgyri field model system. It is used to plot the 
    % dynamics of the Gyorgyri field model in the "Gyorgyri field model_dynamics" function.
    % Input : time (t), concentation (x), rate constants (k1, k2, k3, k4, k5, k6, k7),
    % constant concentrations (A and H)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    xdot = [-k1*H*x(1)*x(2) + k2*H^2*A*x(2) - 2*k3*x(1)^2 - (1/2)*k4*A^(1/2)*H^(1.5)*(C-x(3))*x(1)^(1/2) + (1/2)*k5*x(1)*x(3) ; -k1*H*x(1)*x(2) - k2*A*H^2*x(2) + a*k6*x(3)*x(4) ; k4*A^(0.5)*H^(1.5)*(C-x(3))*x(1)^(0.5) - k5*x(1)*x(3) - a*k6*x(4)*x(3) - b*k7*M*x(3) ; 2*k1*H*x(1)*x(2) + k2*A*H^2*x(2) + k3*x(1)^2 - a*k6*x(4)*x(3)];

end