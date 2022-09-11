function xdot = BZ_reduced_system(t, x, a, b)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function represents the Brusselator system. It is used in other
% functions of this folder to plot the dynamics for example.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    xdot(1,1) = a - (b+1)*x(1) + x(1)^2*x(2);
    xdot(2,1) = b*x(1) - x(1)^2*x(2);
    
%     k1 = 5;
%     
%     xdot(1,1) = a - (b+1)*(x(1)-k1*x(1)) + (x(1)-k1*x(1))^2*(x(2)-k1*x(2));
%     xdot(2,1) = b*(x(1)-k1*x(1)) - (x(1)-k1*x(1))^2*(x(2)-k1*x(2));
end