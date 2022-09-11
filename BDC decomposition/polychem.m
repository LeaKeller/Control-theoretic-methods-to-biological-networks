function X = polychem(B,C,n,m,X,max_itera)

for i = 1:m
    Phi(:,:,i) = eye(n) + B(:,i)*C(i,:);
end

index = 0;
esci = 0;

while esci == 0
    % Range computation
    Y = X;
    for i = 1:m
        Y = [Y Phi(:,:,i)*X];
    end

    % Convex hull
    M=Y;
    K = convhulln(M');
    M=M(:,unique(K(:)));

    % Check equality
        dim_old=size(X);
        dim_old=dim_old(1,2);
        dim_new=size(M);
        dim_new=dim_new(1,2);
        if  dim_new == dim_old
            cambiata = 0;
            for i=1:dim_new
                trovato = 0;
                for j=1:dim_new
                    if(M(:,i)==X(:,j))
                        trovato = 1;
                    end
                end
                if (trovato == 0)
                    cambiata = 1;
                end
            end
            if (cambiata == 0)
                esci = 1;
                fprintf(1, 'M has not changed!\n\n' );
            end
        end
        
        index = index+1;
        if index > max_itera
            esci = 1;
            fprintf(1, 'Maximum number of iterations reached with M always different!\n\n');
        end
        X=M;
end
end