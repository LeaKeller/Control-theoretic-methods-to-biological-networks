%% Functions for Polyhedral procedure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Those functions have been taken from the PhD thesis of Giulia Giordano in
% order to compute BDC decomposition and the matrix of the vertices of the
% polyhedral of the Lyapunov polyhedral function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function X = ini_one(n)
X = [eye(n) -eye(n)];
end