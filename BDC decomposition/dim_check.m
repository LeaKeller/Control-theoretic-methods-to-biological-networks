function [n,m] = dim_check(B,C)

[nB, mB] = size(B);
[mC, nC] = size(C);

if (nC == nB) && (mC == mB)
    n = nB;
    m = mB;
else
    n = 0;
    m = 0;
    fprintf(1, 'ERROR: the input matrices dimension is not consistent!\n\n' );
end
end