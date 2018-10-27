function [U, X, positions, values, tus, err] = g_dla_NoXUpdate(Data, m, positions, values, X, K)

tic;
[n, ~] = size(Data);
err = [];

% number of iterations
% K = 5;

if (sum(positions(:)) == 0)
    % the initialization
    positions = zeros(2, m);
    values = zeros(4, m);

    Z = Data*X';
    scores_nuclear = zeros(n);
    for i = 1:n
        for j = i+1:n
            T = Z([i j], [i j]);
            c1 = norm(T, 'fro')^2/2; c1_2 = c1^2; c2_2 = det(T)^2;
            scores_nuclear(i, j) = sqrt(c1 + sqrt(c1_2 - c2_2)) + sqrt(c1 - sqrt(c1_2 - c2_2)) - trace(T);
        end
    end

    workingX = X;
    for kk = 1:m
        [~, index_nuc] = max(scores_nuclear(:));
        [i_nuc, j_nuc] = ind2sub([n n], index_nuc);

        [Uu, ~, Vv] = svd(Data([i_nuc j_nuc], :)*workingX([i_nuc j_nuc], :)');
        GG = Uu*Vv';

        positions(1, kk) = i_nuc;
        positions(2, kk) = j_nuc;
        values(:, kk) = vec(GG);

        workingX = applyGTransformOnLeft(workingX, i_nuc, j_nuc, values(:, kk));
        Z = Data*workingX';

        for i = [i_nuc j_nuc]
            for j = i+1:n
                T = Z([i j], [i j]);
                c1 = norm(T, 'fro')^2/2; c1_2 = c1^2; c2_2 = det(T)^2;
                scores_nuclear(i, j) = sqrt(c1 + sqrt(c1_2 - c2_2)) + sqrt(c1 - sqrt(c1_2 - c2_2)) - trace(T);
            end
        end

        for j = [i_nuc j_nuc]
            for i = 1:j-1
                T = Z([i j], [i j]);
                c1 = norm(T, 'fro')^2/2; c1_2 = c1^2; c2_2 = det(T)^2;
                scores_nuclear(i, j) = sqrt(c1 + sqrt(c1_2 - c2_2)) + sqrt(c1 - sqrt(c1_2 - c2_2)) - trace(T);
            end
        end
    end
end

err = zeros(K, 1);
for k = 1:K
    the_Data = Data;
    for h = m:-1:1
        the_Data = applyGTransformOnLeftTransp(the_Data, positions(1, h), positions(2, h), values(:, h));
    end
    the_X = X;
    
    for kk = 1:m
        the_Data = applyGTransformOnLeft(the_Data, positions(1, kk), positions(2, kk), values(:, kk));
        
        Z = the_Data*the_X';
        if (kk == 1)
            scores_nuclear = zeros(n);
            for i = 1:n
                for j = i+1:n
                    T = Z([i j], [i j]);
                    c1 = norm(T, 'fro')^2/2; c1_2 = c1^2; c2_2 = det(T)^2;
                    scores_nuclear(i, j) = sqrt(c1 + sqrt(c1_2 - c2_2)) + sqrt(c1 - sqrt(c1_2 - c2_2)) - trace(T);
                end
            end
        else
            for i = [i_nuc j_nuc positions(1, kk) positions(2, kk)]
                for j = i+1:n
                    T = Z([i j], [i j]);
                    c1 = norm(T, 'fro')^2/2; c1_2 = c1^2; c2_2 = det(T)^2;
                    scores_nuclear(i, j) = sqrt(c1 + sqrt(c1_2 - c2_2)) + sqrt(c1 - sqrt(c1_2 - c2_2)) - trace(T);
                end
            end
            
            for j = [i_nuc j_nuc positions(1, kk) positions(2, kk)]
                for i = 1:j-1
                    T = Z([i j], [i j]);
                    c1 = norm(T, 'fro')^2/2; c1_2 = c1^2; c2_2 = det(T)^2;
                    scores_nuclear(i, j) = sqrt(c1 + sqrt(c1_2 - c2_2)) + sqrt(c1 - sqrt(c1_2 - c2_2)) - trace(T);
                end
            end
        end
        
        [~, index_nuc] = max(scores_nuclear(:));
        [i_nuc, j_nuc] = ind2sub([n n], index_nuc);
        
        [Uu, ~, Vv] = svd(the_Data([i_nuc j_nuc], :)*the_X([i_nuc j_nuc], :)');
        GG = Uu*Vv';
        
        positions(1, kk) = i_nuc;
        positions(2, kk) = j_nuc;
        values(:, kk) = vec(GG);
        
        the_X = applyGTransformOnLeft(the_X, positions(1, kk), positions(2, kk), values(:, kk));
    end

    UX = X;
    for h = 1:m
        UX = applyGTransformOnLeft(UX, positions(1, h), positions(2, h), values(:, h));
    end
    err(k) = norm(Data-UX, 'fro')^2/norm(Data, 'fro')^2*100;
end

% explicit dictionary
U = eye(n);
for h = 1:m
    U = applyGTransformOnLeft(U, positions(1, h), positions(2, h), values(:, h));
end
tus = toc;
