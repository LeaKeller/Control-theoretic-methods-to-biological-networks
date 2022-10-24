function nullclines = Nulclines_Oregonator()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the nullclines for the Oregonator system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    A = 0.06;
    B = 0.02;

    k1 = 1.28;
    k2 = 2.4*10^6;
    k3 = 33.6;
    k4 = 2400;
    k5 = 1;

    f= 1;
    q = 1;%(2*k1*k4)/(k2*k3);%1;
    eps = k5/(k3*A);

    % create the mesh for the vector field
    mesh_x = -5:0.5:5;
    mesh_z = -5:0.5:5;
    [x, z] = meshgrid(mesh_x, mesh_z);

    dx = (1/eps).*(((f*z)./(q+x).*(q+x))+x.*(1-x));%(1/eps)*((1-x).*(((f*z)./(1+x))+x));%(1/eps)*((1-x).*(((f*z)./(1+x))+x));
    dz = x-z;
    null_x = x./x;
    null_z = x;

    % cretate the mesh for the streamlines
    [start_x,start_z] = meshgrid(0.5,[0.5,1.5]);

    verts = stream2(x,z,dx,dz,start_x,start_z);

    figure;
    quiver(x,z,dx,dz); hold on;
    plot(mesh_x,mesh_z,'k'); hold on;
    streamline(verts); hold on;
    xline(1,'--r',{'x=1'});
    yline(1,'--r',{'z=1'});
    xlim([0 5.5]);
    ylim([0 5.5]);
    xlabel('x');
    ylabel('z');
end