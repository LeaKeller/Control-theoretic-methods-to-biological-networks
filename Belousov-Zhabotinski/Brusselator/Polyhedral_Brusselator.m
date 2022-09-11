function P = Polyhedral_Brusselator()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    % Boundedness study using Polyhedral function. The code from Guilia PhD
    % thesis is used. Somme text examples from Guilia PhD thesis are also
    % reproduced.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Manually computed B and C 
    
    B = [-1 1 1 0;0 -1 -1 1];
    C = [1 0;1 0;0 1;1 0];

    X1 = polychemtest(B,C,100,1);

    % Using code from Guilia PhD thesis
    
    S = [-1 1 0 ; 0 -1 1];
    g = [1 0 ; 1 2 ; 1 0]; 

    M = SSIM_general(S,g);

    [B, C] = generateBCmatrices(S,g);

    X2 = polychemtest(B,C,100,1);
    
    %% Draw polyhedral
    
    p1 = [2.5353 -2.5353 1.2677 -1.2677];
    p2 = [0 0 1.2677 -1.2677];
    
    poly1 = convhull(p1, p2);
    poly2 = polyshape(p1, p2);
    figure;
    plot(poly2); hold on;
    pts = scatter(p1,p2);
    xlim([-3, 3]);
    ylim([-1.5, 1.5]);

    %% Example 6.4

    x = [1 0 0 -1 0 0 1 -1 -1 1 1 -1];
    y = [0 1 0 0 -1 0 0 0 1 -1 -1 1];
    z = [0 0 1 0 0 -1 1 -1 0 0 1 -1];
    figure;
    pts = scatter3(x,y,z);
    figure;
    [k,av] = convhull(x,y,z);
    trisurf(k,x,y,z,'FaceColor','cyan');
    figure;
    plot3(x,y,z);

    %% Example 6.3.2 Giulia Giordano 

    X = [1.0000 0 -1.0000 0 0.5000 -0.5000 0.0000 -0.0000 1.0000 -1.0000];
    Y = [0 0 0 0 0 0 1.0000 -1.0000 -2.0000 2.0000];
    Z = [0 1.0000 0 -1.0000 1.0000 -1.0000 1.0000 -1.0000 0 0];
    figure;
    points = scatter3(X,Y,Z);
    figure;
    [k1,av1] = convhull(X,Y,Z);
    trisurf(k1,X,Y,Z,'FaceColor','cyan');
    figure;
    plot3(X,Y,Z);
end