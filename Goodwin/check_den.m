%% Functions for steady state influcence matrix

function outcome = check_den(B,C,epsitol)

[n m] = size(B);

% Flags
flag_p = 0; % positive
flag_n = 0; % negative
flag_z = 0; % zero

% Compute det(-BDC) for all the possible choices
% of the diagonal entries Dk of the diagonal matrix D
% with each Dk being either 0 or 1
for i = 0: 2^m - 1

    % Create the binary vector
    c = double(bitget(uint64(i), 1:1:m));
    % Build the diagonal matrix
    D = diag(c);
    % Compute det(-J), J=BDC
    den = det(-B*D*C);

    if abs(den) < epsitol % The determinant is zero
        flag_z = flag_z + 1;
    elseif den < 0 % The determinant is negative
        flag_n = flag_n + 1;
        if flag_p > 0
            break
        end
    elseif den > 0 % The determinant is positive
        flag_p = flag_p + 1;
        if flag_n > 0
            break
        end
    end
end

if (flag_n == 0 && flag_p == 0)
    % The determinant is structurally zero
    outcome = 0;
elseif (flag_p == 0 && det(-B*C) < 0)
    % The determinant is structurally negative
    outcome = -1;
elseif (flag_n == 0 && det(-B*C) > 0)
    % The determinant is structurally positive
    outcome = 1;
else
    % The sign is indeterminate
    outcome = -99;
end

end