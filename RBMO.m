%% from ref.RBMO
function [BestValue, Xfood, Conv, result] = RBMO(N, T, Xmin, Xmax, D, fobj)
    Xfood = zeros(1, D);
    BestValue = inf;
    Conv = zeros(1, T);
    fitness = inf(N, 1);
    Epsilon = 0.5;
    
    X = initialization(N, D, Xmax, Xmin);
    
    X_old = X;
    for i = 1:N
        fitness_old(i, 1) = fobj(X_old(i, :));
    end
  
    t = 1; 
    
    while t < T+1
        t

        for i = 1:N

            p = randi([2, 5]);
            selected_index_p = randperm(N, p);
            Xp = X(selected_index_p, :);
            Xpmean = mean(Xp);

            q = randi([10, N]);
            selected_index_q = randperm(N, q);
            Xq = X(selected_index_q, :);
            Xqmean = mean(Xq);

            A = randperm(N);
            R1 = A(1);
            if rand < Epsilon
                X(i, :) = X(i, :) + (Xpmean - X(R1, :)) .* rand; % Eq. (3)
            else
                X(i, :) = X(i, :) + (Xqmean - X(R1, :)) .* rand; % Eq. (4)
            end
        end
        

        X = boundaryCheck(X, Xmin, Xmax);
        
        for i = 1:N
            fitness(i, 1) = fobj(X(i, :));
            if fitness(i, 1) < BestValue
                BestValue = fitness(i, 1);
                Xfood = X(i, :);
            end
        end
        

        [fitness, X, fitness_old, X_old] = Food_storage(fitness, X, fitness_old, X_old); % Eq. (7)
        
        CF = (1 - t / T)^(2 * t / T);
        
        for i = 1:N
            p = randi([2, 5]);
            selected_index_p = randperm(N, p);
            Xp = X(selected_index_p, :);
            Xpmean = mean(Xp);

            q = randi([10, N]);
            selected_index_q = randperm(N, q);
            Xq = X(selected_index_q, :);
            Xqmean = mean(Xq);

            if rand() < Epsilon
                X(i, :) = Xfood + CF * (Xpmean - X(i, :)) .* randn(1, D); % Eq. (5)
            else
                X(i, :) = Xfood + CF * (Xqmean - X(i, :)) .* randn(1, D); % Eq. (6)
            end
        end

        X = boundaryCheck(X, Xmin, Xmax);

        for i = 1:N
            fitness(i, 1) = fobj(X(i, :));
            if fitness(i, 1) < BestValue
                BestValue = fitness(i, 1);
                Xfood = X(i, :);
            end
        end

        [fitness, X, fitness_old, X_old] = Food_storage(fitness, X, fitness_old, X_old); % Eq. (7)
        
        Conv(t) = BestValue;
        result(t, :)= Xfood;
        t = t + 1;
    end
end


function [fit, X, fit_old, X_old] = Food_storage(fit, X, fit_old, X_old)
    Inx = (fit_old < fit);
    Indx = repmat(Inx, 1, size(X, 2));
    X = Indx .* X_old + ~Indx .* X;
    fit = Inx .* fit_old + ~Inx .* fit;
    fit_old = fit;
    X_old = X;
end


function X = boundaryCheck(X, Xmin, Xmax)
    for i = 1:size(X, 1)
        FU = X(i, :) > Xmax;
        FL = X(i, :) < Xmin;
        X(i, :) = (X(i, :) .* (~(FU + FL))) + Xmax .* FU + Xmin .* FL;
    end
end