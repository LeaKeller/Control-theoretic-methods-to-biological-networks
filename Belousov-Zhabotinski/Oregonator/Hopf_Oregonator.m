function H = Hopf_Oregonator()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the different eigenvalues for the Oregonator
% system. First, stability is checked by plotting the eigenvalues and the
% unit circle. Second, the Hopf bifurcations are performed regardinf the
% different equilibria.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    syms x z
    
    A = 0.06;
    B = 0.02;
    
    k1 = 1.28;
    k2 = 2.4*10^6;
    k3 = 33.6;
    k4 = 2400;
    k5 = 1;
    
    q = (2*k1*k4)/(k2*k3);
    e = k5/(k3*A);
    f = 1;
    
    f1 = (q/e)*((f*z)/(q+x)) - (1/e)*x*((f*z)/(q+x)) + (1/e)*x*(1-x);
    f2 = x-z;
    
    sol = vpasolve([f1 == 0 f2 == 0], [x,z]); % equilibrium points
    
    Df1x = diff(f1,x);
    Df1z = diff(f1,z);
    
    Df1x0 = [];
    Df1z0 = [];

    % Compute the equilibria using Jacobian matrix for the 2 variables reduced
    % model
    
    for i=1:size(sol.x,1)
        Df1x0 = [Df1x0 subs(Df1x, [x z], [sol.x(i), sol.z(i)])];
        Df1z0 = [Df1z0 subs(Df1z, [x z], [sol.x(i), sol.z(i)])];
    end
    
    % Compute the Jacobian matrices
    
    J1 = [Df1x0(1) Df1z0(1) ; 1 -1];
    J2 = [Df1x0(1) Df1z0(2) ; 1 -1];
    J3 = [Df1x0(1) Df1z0(3) ; 1 -1];
    J4 = [Df1x0(2) Df1z0(2) ; 1 -1];
    J5 = [Df1x0(3) Df1z0(3) ; 1 -1];
    J6 = [Df1x0(3) Df1z0(3) ; 1 -1];
    
    % Check that the eigenvalues are in the unit circle
    
    eigenval1 = eig(J1);     
    eigenval2 = eig(J2); 
    eigenval3 = eig(J3); 
    eigenval4 = eig(J4);     
    eigenval5 = eig(J5); 
    eigenval6 = eig(J6);

    l_eig1_real = [real(eigenval1(1)) real(eigenval2(1)) real(eigenval3(1)) real(eigenval4(1)) real(eigenval5(1)) real(eigenval6(1))];
    l_eig2_real = [real(eigenval1(2)) real(eigenval2(2)) real(eigenval3(2)) real(eigenval4(2)) real(eigenval5(2)) real(eigenval6(2))];
    l_eig1_imag = [imag(eigenval1(1)) imag(eigenval2(1)) imag(eigenval3(1)) imag(eigenval4(1)) imag(eigenval5(1)) imag(eigenval6(1))];
    l_eig2_imag = [imag(eigenval1(2)) imag(eigenval2(2)) imag(eigenval3(2)) imag(eigenval4(2)) imag(eigenval5(2)) imag(eigenval6(2))];
    
    figure;
    for i=1:size(l_eig1_real,1)*size(l_eig1_real,2)
        scatter([l_eig1_real(i), l_eig2_real(i)], [l_eig1_imag(i), l_eig2_imag(i)]); hold on; 
    end
    circle = @(x,y) x.^2 + y.^2 -1;
    ezplot(circle); hold on; % unit circle
    axis([-3 3 -3 3]);
    
    %Hopf bifurcation
    
    TrJ = (2*f*z-(1-2*x)*(q+x))/(e*(q+x));
    DetJ = -((-2*f*z+(1-2*x)*(q+x))/(e*(q+x))+(f*(q-x))/(e*(q+x)));
    
    trJ = [];
    detJ =[];

    trJ1 = [trJ subs(TrJ, [x z], [sol.x(1), sol.z(1)])]; % substitute the value for TrJ at equilibria
    trJ2 = [trJ subs(TrJ, [x z], [sol.x(1), sol.z(2)])];
    trJ3 = [trJ subs(TrJ, [x z], [sol.x(1), sol.z(3)])];
    trJ4 = [trJ subs(TrJ, [x z], [sol.x(2), sol.z(2)])];
    trJ5 = [trJ subs(TrJ, [x z], [sol.x(3), sol.z(2)])];
    trJ6 = [trJ subs(TrJ, [x z], [sol.x(3), sol.z(3)])];
    
    trJ = [trJ1 trJ2 trJ3 trJ4 trJ5 trJ6];
  
    detJ1 = [detJ subs(DetJ, [x z], [sol.x(1), sol.z(1)])]; % substitute the value for DetJ at equilibria
    detJ2 = [detJ subs(DetJ, [x z], [sol.x(1), sol.z(2)])];
    detJ3 = [detJ subs(DetJ, [x z], [sol.x(1), sol.z(3)])];
    detJ4 = [detJ subs(DetJ, [x z], [sol.x(2), sol.z(2)])];
    detJ5 = [detJ subs(DetJ, [x z], [sol.x(3), sol.z(2)])];
    detJ6 = [detJ subs(DetJ, [x z], [sol.x(3), sol.z(3)])];
    
    detJ = [detJ1 detJ2 detJ3 detJ4 detJ5 detJ6];
    
    % Hopf bifurcation with e as bifurcation parameter
    
    syms e q f
    
    f1 = (q/e)*((f*z)/(q+x)) - (1/e)*x*((f*z)/(q+x)) + (1/e)*x*(1-x);
    f2 = x-z;

    sol = vpasolve([f1 == 0 f2 == 0], [x,z]); % equilibrium points

    Df1x = diff(f1,x);
    Df1z = diff(f1,z);
    
    Df1x_1 = subs(Df1x, [x z], [sol.x(1), sol.z(1)]);
    Df1z_1 = subs(Df1z, [x z], [sol.x(1), sol.z(2)]);
    Df1x_2 = subs(Df1x, [x z], [sol.x(1), sol.z(3)]);
    Df1z_2 = subs(Df1z, [x z], [sol.x(2), sol.z(2)]);
    Df1x_3 = subs(Df1x, [x z], [sol.x(3), sol.z(2)]);
    Df1z_3 = subs(Df1z, [x z], [sol.x(3), sol.z(3)]);
    
    L1_1 = -(1/2)*((Df1x_1-1)+4*sqrt((1-Df1x_1)^2-4*(Df1z_1)));
    L2_1= -(1/2)*((Df1x_1-1)-4*sqrt((1-Df1x_1)^2-4*(Df1z_1)));
    L1_2 = -(1/2)*((Df1x_2-1)+4*sqrt((1-Df1x_2)^2-4*(Df1z_2)));
    L2_2 = -(1/2)*((Df1x_2-1)-4*sqrt((1-Df1x_2)^2-4*(Df1z_2)));
    L1_3 = -(1/2)*((Df1x_3-1)+4*sqrt((1-Df1x_3)^2-4*(Df1z_3)));
    L2_3= -(1/2)*((Df1x_3-1)-4*sqrt((1-Df1x_3)^2-4*(Df1z_3)));
    
    n = 14;
    f=1; 
    q = (2*(f+1)+(2*f+2)^2-4*(f^2+6*f-3))/2;
    e = [];
    for j = -1:n-2
        e = [e, 1.5*j];
    end

    c = ['r','g','b','c','m','y','k','w',[0.6350 0.0780 0.1840],[0.3010 0.7450 0.9330],[0.4660 0.6740 0.1880],[0.4940 0.1840 0.5560],[0.9290 0.6940 0.1250],[0.8500 0.3250 0.0980],[0 0.4470 0.7410]];
    
    figure;
    for i = 1:n
        txt = ['e =', num2str(e(i))];
        l1_1(i) = 1/2 - 2*((1/e(i) - 1)^2 - (4*f)/e(i))^(1/2) - 1/(2*e(i));
        l2_1(i) = (f*(0.5*f + 0.5*q + 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) - 0.5))/(2*e(i)*(0.5*f - 0.5*q + 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) - 0.5)) - (0.5*f + 0.5*q + 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) - 0.5)/(2*e(i)) - (0.5*f + 0.5*q + 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) + 0.5)/(2*e(i)) - 2*(((0.5*f + 0.5*q + 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) - 0.5)/e(i) + (0.5*f + 0.5*q + 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) + 0.5)/e(i) - (f*(0.5*f + 0.5*q + 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) - 0.5))/(e(i)*(0.5*f - 0.5*q + 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) - 0.5)) + (f*(0.5*f + 0.5*q + 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) - 0.5)^2)/(e(i)*(0.5*f - 0.5*q + 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) - 0.5)^2) + (f*q*(0.5*f + 0.5*q + 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) - 0.5))/(e(i)*(0.5*f - 0.5*q + 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) - 0.5)^2) - 1)^2 + (4*f*(0.5*f + 0.5*q + 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) - 0.5))/(e(i)*(0.5*f - 0.5*q + 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) - 0.5)) + (4*f*q)/(e(i)*(0.5*f - 0.5*q + 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) - 0.5)))^(1/2) - (f*(0.5*f + 0.5*q + 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) - 0.5)^2)/(2*e(i)*(0.5*f - 0.5*q + 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) - 0.5)^2) - (f*q*(0.5*f + 0.5*q + 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) - 0.5))/(2*e(i)*(0.5*f - 0.5*q + 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) - 0.5)^2) + 1/2;
        scatter([real(l1_1(i)),real(l2_1(i))], [imag(l1_1(i)),imag(l2_1(i))], c(i), 'DisplayName', txt); hold on;
    end
    title('f=1, q=2');
    plot(real(l1_1),imag(l1_1),'DisplayName','\lambda_1');
    plot(real(l2_1),imag(l2_1),'DisplayName','\lambda_2');
    xlabel('\Re(\lambda)');
    ylabel('\Im(\lambda)');
    legend show;
    
    figure;
    for i = 1:n
        txt = ['e =', num2str(e(i))];
        l1_2(i) = (f*(0.5*f + 0.5*q + 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) - 0.5))/(2*e(i)*(0.5*f - 0.5*q + 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) - 0.5)) - (0.5*f + 0.5*q + 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) - 0.5)/(2*e(i)) - (0.5*f + 0.5*q + 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) + 0.5)/(2*e(i)) - 2*(((0.5*f + 0.5*q + 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) - 0.5)/e(i) + (0.5*f + 0.5*q + 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) + 0.5)/e(i) - (f*(0.5*f + 0.5*q + 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) - 0.5))/(e(i)*(0.5*f - 0.5*q + 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) - 0.5)) + (f*(0.5*f + 0.5*q + 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) - 0.5)^2)/(e(i)*(0.5*f - 0.5*q + 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) - 0.5)^2) + (f*q*(0.5*f + 0.5*q + 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) - 0.5))/(e(i)*(0.5*f - 0.5*q + 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) - 0.5)^2) - 1)^2 + (4*f*(0.5*f + 0.5*q + 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) - 0.5))/(e(i)*(0.5*f - 0.5*q + 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) - 0.5)) + (4*f*q)/(e(i)*(0.5*f - 0.5*q + 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) - 0.5)))^(1/2) - (f*(0.5*f + 0.5*q + 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) - 0.5)^2)/(2*e(i)*(0.5*f - 0.5*q + 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) - 0.5)^2) - (f*q*(0.5*f + 0.5*q + 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) - 0.5))/(2*e(i)*(0.5*f - 0.5*q + 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) - 0.5)^2) + 1/2;
        l2_2(i) = 2*(((0.5*f + 0.5*q + 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) - 0.5)/e(i) + (0.5*f + 0.5*q + 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) + 0.5)/e(i) - (f*(0.5*f + 0.5*q + 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) - 0.5))/(e(i)*(0.5*f - 0.5*q + 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) - 0.5)) + (f*(0.5*f + 0.5*q + 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) - 0.5)^2)/(e(i)*(0.5*f - 0.5*q + 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) - 0.5)^2) + (f*q*(0.5*f + 0.5*q + 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) - 0.5))/(e(i)*(0.5*f - 0.5*q + 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) - 0.5)^2) - 1)^2 + (4*f*(0.5*f + 0.5*q + 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) - 0.5))/(e(i)*(0.5*f - 0.5*q + 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) - 0.5)) + (4*f*q)/(e(i)*(0.5*f - 0.5*q + 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) - 0.5)))^(1/2) - (0.5*f + 0.5*q + 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) - 0.5)/(2*e(i)) - (0.5*f + 0.5*q + 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) + 0.5)/(2*e(i)) + (f*(0.5*f + 0.5*q + 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) - 0.5))/(2*e(i)*(0.5*f - 0.5*q + 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) - 0.5)) - (f*(0.5*f + 0.5*q + 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) - 0.5)^2)/(2*e(i)*(0.5*f - 0.5*q + 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) - 0.5)^2) - (f*q*(0.5*f + 0.5*q + 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) - 0.5))/(2*e(i)*(0.5*f - 0.5*q + 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) - 0.5)^2) + 1/2;
        scatter([real(l1_2(i)),real(l2_2(i))], [imag(l1_2(i)),imag(l2_2(i))], c(i), 'DisplayName', txt); hold on;
    end
    title('f=1, q=2');
    plot(real(l1_2),imag(l1_2),'DisplayName','\lambda_1');
    plot(real(l2_2),imag(l2_2),'DisplayName','\lambda_2');
    xlabel('\Re(\lambda)');
    ylabel('\Im(\lambda)');
    legend show;
    
    figure;
    for i = 1:n
        txt = ['e =', num2str(e(i))];
        l1_3(i) = 1/2 - (0.5*f + 0.5*q - 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) - 0.5)/(2*e(i)) - (0.5*f + 0.5*q - 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) + 0.5)/(2*e(i)) - (f*(0.5*f + 0.5*q - 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) - 0.5))/(2*e(i)*(0.5*q - 0.5*f + 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) + 0.5)) - (f*(0.5*f + 0.5*q - 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) - 0.5)^2)/(2*e(i)*(0.5*q - 0.5*f + 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) + 0.5)^2) - (f*q*(0.5*f + 0.5*q - 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) - 0.5))/(2*e(i)*(0.5*q - 0.5*f + 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) + 0.5)^2) - 2*(((0.5*f + 0.5*q - 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) - 0.5)/e(i) + (0.5*f + 0.5*q - 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) + 0.5)/e(i) + (f*(0.5*f + 0.5*q - 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) - 0.5))/(e(i)*(0.5*q - 0.5*f + 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) + 0.5)) + (f*(0.5*f + 0.5*q - 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) - 0.5)^2)/(e(i)*(0.5*q - 0.5*f + 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) + 0.5)^2) + (f*q*(0.5*f + 0.5*q - 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) - 0.5))/(e(i)*(0.5*q - 0.5*f + 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) + 0.5)^2) - 1)^2 - (4*f*(0.5*f + 0.5*q - 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) - 0.5))/(e(i)*(0.5*q - 0.5*f + 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) + 0.5)) - (4*f*q)/(e(i)*(0.5*q - 0.5*f + 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) + 0.5)))^(1/2);
        l2_3(i) = 2*(((0.5*f + 0.5*q - 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) - 0.5)/e(i) + (0.5*f + 0.5*q - 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) + 0.5)/e(i) + (f*(0.5*f + 0.5*q - 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) - 0.5))/(e(i)*(0.5*q - 0.5*f + 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) + 0.5)) + (f*(0.5*f + 0.5*q - 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) - 0.5)^2)/(e(i)*(0.5*q - 0.5*f + 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) + 0.5)^2) + (f*q*(0.5*f + 0.5*q - 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) - 0.5))/(e(i)*(0.5*q - 0.5*f + 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) + 0.5)^2) - 1)^2 - (4*f*(0.5*f + 0.5*q - 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) - 0.5))/(e(i)*(0.5*q - 0.5*f + 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) + 0.5)) - (4*f*q)/(e(i)*(0.5*q - 0.5*f + 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) + 0.5)))^(1/2) - (0.5*f + 0.5*q - 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) - 0.5)/(2*e(i)) - (0.5*f + 0.5*q - 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) + 0.5)/(2*e(i)) - (f*(0.5*f + 0.5*q - 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) - 0.5))/(2*e(i)*(0.5*q - 0.5*f + 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) + 0.5)) - (f*(0.5*f + 0.5*q - 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) - 0.5)^2)/(2*e(i)*(0.5*q - 0.5*f + 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) + 0.5)^2) - (f*q*(0.5*f + 0.5*q - 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) - 0.5))/(2*e(i)*(0.5*q - 0.5*f + 0.5*(f^2 + 6.0*f*q - 2.0*f + q^2 + 2.0*q + 1.0)^(1/2) + 0.5)^2) + 1/2;
        scatter([real(l1_3(i)),real(l2_3(i))], [imag(l1_3(i)),imag(l2_3(i))], c(i), 'DisplayName', txt); hold on;
    end
    title('f=1, q=2');
    plot(real(l1_3),imag(l1_3),'DisplayName','\lambda_1');
    plot(real(l2_3),imag(l2_3),'DisplayName','\lambda_2');
    xlabel('\Re(\lambda)');
    ylabel('\Im(\lambda)');
    legend show;

    % compute the eigenvalues of the Jacobian 
    % First compute the characteristic polynomial 
    syms lam
    Jac_l = [b-1-lam a^2 ; -b -a^2-lam];
    lamP = charpoly(Jac,lam);
    eigenvalues = eig(Jac_l);

    format rat

    n = 14;
    a=3; 
    b = [];
    for j = 1:n
        b = [b, 1.5*j];
    end
    c = ['r','g','b','c','m','y','k','w',[0.6350 0.0780 0.1840],[0.3010 0.7450 0.9330],[0.4660 0.6740 0.1880],[0.4940 0.1840 0.5560],[0.9290 0.6940 0.1250],[0.8500 0.3250 0.0980],[0 0.4470 0.7410]];
    
    figure;
    for i = 1:n
        txt = ['b =', num2str(b(i))];
        l1(i) = b(i)/2 - (-(- a^2 + 2*a + b(i) - 1)*(a^2 + 2*a - b(i) + 1))^(1/2)/2 - a^2/2 - 1/2;
        l2(i) = b(i)/2 + (-(- a^2 + 2*a + b(i) - 1)*(a^2 + 2*a - b(i) + 1))^(1/2)/2 - a^2/2 - 1/2;
        scatter([real(l1(i)),real(l2(i))], [imag(l1(i)),imag(l2(i))], c(i), 'DisplayName', txt); hold on;
    end
    title('a=3');
    plot(real(l1),imag(l1),'DisplayName','\lambda_1');
    plot(real(l2),imag(l2),'DisplayName','\lambda_2');
    xlabel('\Re(\lambda)');
    ylabel('\Im(\lambda)');
    legend show;
end