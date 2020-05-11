function [S_k,S_k_idx] = generate_S2(L,R)
%% Function takes 
% L: Number of Lags
% R: Number of Regions
% as arguments first generate a LR^2 times LR^2 symmetrix matrix of binary
% elements to be used for the ICAR prior.
% The construction of the matrix follows the methodology of Chiang Et al.
% The function returns a matrix S_k, which stores the sum of neighbors of
% the kth entry, and structure S_k_idx, which stores the actual indices of
% the neighbors
S_k = zeros(1,L*R^2);
S_k_idx = cell(1,L*R^2);
% of the matrix
block_mat = kron(eye(L*R),ones(R,R));
A = zeros(L*R^2,L*R^2);
if L > 1
    for j = 1:R
        for j_prime = 1:R
            for l = 1:(L-1)
                A((j-1)*L*R+l*R+j_prime, (j-1)*L*R+j_prime) = 1;
            end
        end
    end
else
    for j = 1:R
        for j_prime = 1:R
            A((j-1)*L*R+R+j_prime, (j-1)*L*R+j_prime) = 1;
        end
    end
end
A = A(1:(L*R^2),:);
S = block_mat + A + A';

% Now we convert S into a sparse matrix
S = sparse(S);
for i = 1:(L*R^2)
    S_k_idx{i}= setdiff(find(S(i,:)==1),i);
    S_k(i) = length(S_k_idx{i});
end
end
