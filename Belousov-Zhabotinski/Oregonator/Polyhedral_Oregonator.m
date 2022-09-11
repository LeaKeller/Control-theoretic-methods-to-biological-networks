function P = Polyhedral_Oregonator()
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    % Boundedness study using Polyhedral function. The code from Guilia PhD
    % thesis is used. Somme text examples from Guilia PhD thesis are also
    % reproduced.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    A = 0.06;
    B = 0.02;
    
    k1 = 1.28;
    k2 = 2.4*10^6;
    k3 = 33.6;
    k4 = 2400;
    k5 = 1;
    
    f = 1;
    
    e = k5/k3*A;
    ee = 2*k4*k5/(k2*k3*A);
    q = 2*k1*k4/(k2*k3);

    % Using code from Guilia PhD thesis
    
    S = [1 0 0 0;0 1 1 0;0 0 -1 1];
    g = [1 2 0;1 2 0;0 0 3;1 0 0];

    D = [1 0 0 0 0 0;0 2 0 0 0 0;0 0 3 0 0 0;0 0 0 4 0 0;0 0 0 0 5 0;0 0 0 0 0 6];
    
    M = SSIM_general(S,g);

    [B, C] = generateBCmatrices(S,g);

    X2 = polychemtest(B,C,100,1); 
    
    %% Draw polyhedral
    
    p1 = [2.5353 -2.5353 0 0 0 0];
    p2 = [0 0 2.5353 -2.5353 1.2677 -1.2677];
    p3 = [0 0 0 0 1.2677 -1.2677];
    
    figure;
    poly = polyshape();
    plot3(p1, p2, p3); hold on;
    pts = scatter3(p1,p2,p3);
    grid on;
    figure;
    [k1,av1] = convhull(p1,p2,p3);
    trisurf(k1,p1,p2,p3,'FaceColor','cyan');
end