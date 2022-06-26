function B = Brusselator_dynamics()

    %% Values taken from Yishao Zhou paper
    %%
    
    a = 3;
    b = [0.3, 16.5];
    
    tfinal = 100;
    
    x01 = 0.1;
    x02 = 0.1;
    
    for i = 1:length(b)
        %x_eq = [a; b(i)/a];
        %x_0 = x_eq + 0.01;
        [t, x] = ode45(@(t, x) BZ_reduced_system(t, x, a, b(i)) , [0 tfinal], [x01 ; x02]); %x_0
    
    x1 = x(:,1);
    x2 = x(:,2);
    
    figure;
    plot(x1,x2); hold on;
    %scatter(x_0(1),x_0(2));
    %scatter(x_eq(1),x_eq(2));
    xlabel('x_1');
    ylabel('x_2');
    grid on;
    title(['a=3, ','b=', num2str(b(i))]);
    
    figure;
    plot(t,x1); hold on; 
    plot(t,x2);
    
    xlabel('t');
    legend('x_1','x_2');
    grid on; 
    title(['a=3, ','b=', num2str(b(i))]);
    end
end