function M = SSIM_general(S,g,comb)

% Compute matrices B and C of the BDC-decomposition
if(exist('comb', 'var'))
    [B,C] = generateBCmatrices(S,g,comb);
else
    [B,C] = generateBCmatrices(S,g);
end

% Compute the structural influence matrix
M = SignMatrixGeneration(B,C);

end