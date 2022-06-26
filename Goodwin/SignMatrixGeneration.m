%% Functions for steady state influcence matrix

function signmatrix = SignMatrixGeneration(B,C)

% Numerical tolerance
epsitol = 10^(-10);

% Influence matrix initialisation
n = size(B,1);
signmatrix = zeros(n,n);

% Check the sign of det(-J), J=BDC
cd = check_den(B,C,epsitol);
if cd ~= -99
    if cd == 1
        fprintf(1, '\n det(-J) is structurally positive!\n Thus, if the system is bounded, the equilibrium is unique!\n');
    elseif cd == -1
        fprintf(1, '\n WARNING: det(-J) is structurally negative!\n');
    elseif cd == 0
        fprintf(1, '\n WARNING: det(-J) is structurally zero!\n');
    end
else
    fprintf(1, '\n WARNING: det(-J) is not sign definite!\n');
end

% Fill in the influence matrix,
% namely check the sign of det[-J -E_j; H_i 0], J=BDC,
% for all possible combinations, where:
% E_j has a 1 in the j-th position and zero elsewhere
% H_i has a 1 in the i-th position and zero elsewhere
for j = 1:n
    E = zeros(n,1);
    E(j,1) = 1;
    
    for i = 1:n
        H = zeros(1,n);
        H(1,i) = 1;
                
        cn = check_num(B,C,H,E,epsitol);
        if cn ~= -99
            if cn > 0
                signmatrix(i,j)=1;
            elseif cn < 0
                signmatrix(i,j)= -1;
            end
        else
            signmatrix(i,j)= 2;
        end
        
    end
end

fprintf(1, '\n Here is the influence matrix:\n');
signmatrix

end