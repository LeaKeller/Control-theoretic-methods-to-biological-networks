%% Functions for steady state influcence matrix

function [B,C] = generateBCmatrices(S,g,comb)

B=[];
C=[];
n = size(S,1);

if(~exist('comb', 'var'))
    for i=1:size(S,2)
        [deriv_i, deriv_j]=find(g(i,:)~=0);
        for j=1:size(deriv_i,2)
            B=[B S(:,i)*sign(g(i,deriv_j(1,j)))];
            C=[C; zeros(1,n)];
            C(end,abs(g(i,deriv_j(1,j))))=1;
        end
    end
else
    for i=1:size(S,2)
        [deriv_i, deriv_j]=find(g(i,:)~=0);
        if comb(i)~=0
            B=[B S(:,i)];
            C=[C; zeros(1,n)];
            for j=1:size(deriv_i,2)
                C(end,abs(g(i,deriv_j(1,j))))=1*sign(g(i,deriv_j(1,j)));
            end
        else
            for j=1:size(deriv_i,2)
                B=[B S(:,i)*sign(g(i,deriv_j(1,j)))];
                C=[C; zeros(1,n)];
                C(end,abs(g(i,deriv_j(1,j))))=1;
            end
        end
    end
end

end