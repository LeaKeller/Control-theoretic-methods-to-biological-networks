%% 3 variables model

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the function used in order to plot the the dynamics for the
% Goodwin model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ydot = system(t, y, k1, k2, k3, k4, k5, k6, ki, n) 
 
  ydot(1,1) = k1*(ki^n./(ki^n+y(3)^n)) - k2*y(1);
  ydot(2,1) = k3*y(1) - k4*y(2);
  ydot(3,1) = k5*y(2) - k6*y(3);
 
end