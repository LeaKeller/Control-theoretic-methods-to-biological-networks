function Xdot = Goodwin(t,X,k1,k2,k3,k4,k5,k6,ki,n) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function represents the Goodwin system. It is used in other
% functions of this folder to for instance plot the dynamics.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Xdot = [k1*(ki^n/(ki^n+X(3)^n)) - k2*X(1); k3*X(1) - k4*X(2); k5*X(2)- k6*X(3)]; % non reduced system
%     alpha = k2/(k1*k3*(k5/ki))^(1/3);
%     beta = k4/(k1*k3*(k5/ki))^(1/3);
%     gamma = k6/(k1*k3*(k5/ki))^(1/3);
%     Xdot = [1/(1+X(3)^n - alpha*X(1)) ; X(1) - beta*X(2) ; X(2) - gamma*X(3)]; % reduced system
end