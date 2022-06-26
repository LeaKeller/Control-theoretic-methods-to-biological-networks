function xdot = BZ_reduced_system(t, x, a, b)
    xdot(1,1) = a - (b+1)*x(1) + x(1)^2*x(2);
    xdot(2,1) = b*x(1) - x(1)^2*x(2);
    
%     k1 = 5;
%     
%     xdot(1,1) = a - (b+1)*(x(1)-k1*x(1)) + (x(1)-k1*x(1))^2*(x(2)-k1*x(2));
%     xdot(2,1) = b*(x(1)-k1*x(1)) - (x(1)-k1*x(1))^2*(x(2)-k1*x(2));
end