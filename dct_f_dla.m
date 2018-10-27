function [U, S, supportS, X, tus, err] = dct_f_dla(Data, k0, m, s)

tic;
[n, ~] = size(Data);
U = dctmtx(n);

S = zeros(n,m);
totalsupport = [];
for j = 1:round(s*m/n)+1
    totalsupport = [totalsupport 1:n];
end
totalsupport = totalsupport(randsample(length(totalsupport), length(totalsupport)));

supportS = zeros(s, m);
for i = 1:m
    support = totalsupport(1:s);
    totalsupport(1:s) = [];
    support = sort(support);
    support = unique(support);
    
    index = 0;
    while (length(support) < s)
        index = index + 1;
        support = [support totalsupport(index)];
        support = unique(support);
    end
    
    supportS(:, i) = support;
    S(support,i) = randn(s, 1);
    S(:, i) = S(:,i)/norm(S(:,i));
end
X = omp(S'*U'*Data, S'*S, k0);

for i = 1:m
    notsupport = setdiff(1:m, i);
    J = find(X(i, :));

    R = U'*Data(:, J) - S(:, notsupport)*X(notsupport, J);
    support = supportS(:, i);
    [u, sigma, v] = svds(R(support, :), 1);

    S(support, i) = u;
    X(i, J) = sigma*v';
end

% number of iterations
K = 200;

err = zeros(K, 1);
for k = 1:K
    for i = 1:m
        notsupport = [1:i-1 i+1:m];
        J = find(X(i, :));

        R = U'*Data(:, J) - S(:, notsupport)*X(notsupport, J);
        support = supportS(:, i);
        [u, sigma, v] = svds(R(support, :), 1);

        S(support, i) = u;
        X(i, J) = sigma*v';
    end
    
    X = omp(S'*U'*Data, S'*S, k0); 
    
    err(k) = norm(Data-U*S*X, 'fro')^2/norm(Data, 'fro')^2*100;
end
tus = toc;
