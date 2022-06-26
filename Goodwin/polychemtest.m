function X = polychemtest(B,C,max_itera,dual_procedure)
[n,m] = dim_check(B,C);
if dual_procedure == 1
    [B,C] = dual(B,C);
end
X = ini_one(n);
X = polychem(B,C,n,m,X,max_itera);
end