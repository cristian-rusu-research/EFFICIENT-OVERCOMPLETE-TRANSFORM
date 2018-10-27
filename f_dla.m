function [U, S, X, positions, values, tus, err] = f_dla(Data, k0, g, m, U, S, supportS, X)

tic;
[n, ~] = size(Data);

% number of iterations
K = 200;

positions = zeros(2, g);
values = zeros(4, g);

for i = 1:1
    [U, ~, positions, values, ~, errg] = g_dla_NoXUpdate(Data, g, positions, values, S*X, 3);
    X = omp(S'*U'*Data, S'*S, k0);
end

err = zeros(K, 1);
for k = 1:K
    for kk = 1:3
        for i = 1:m
            notsupport = [1:i-1 i+1:m];
            J = find(X(i, :));

            R = U'*Data(:, J) - S(:, notsupport)*X(notsupport, J);
            support = supportS(:, i);
            [u, sigma, v] = svds(R(support, :), 1);

            S(support, i) = u;
            X(i, J) = sigma*v';
        end

        [U, ~, positions, values, ~, errg] = g_dla_NoXUpdate(Data, g, positions, values, S*X, 1);
    end
    
    X = omp(S'*U'*Data, S'*S, k0);
    
    err(k) = norm(Data-U*S*X, 'fro')^2/norm(Data, 'fro')^2*100;
end
tus = toc;
