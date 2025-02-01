function [X, Y] = generate_nonlinear_system(n, seed, mode)
    rng(seed);

    n = n + 100;

    if mode == "l-uni"
        % n = 20
        X = zeros(n, 1);
        Y = zeros(n, 1);
        X(1) = rand();
        Y(1) = rand();
        gamma_yx = 0.32;
        r_x = 3.7; r_y = 3.8;
        for i = 2:n
            X(i) = X(i - 1) * (r_x - r_x * X(i - 1));
            Y(i) = Y(i - 1) * (r_y - r_y * Y(i - 1) - gamma_yx * X(i - 1));
        end

    elseif mode == "l-bi"
        % n = 20
        X = zeros(n, 1);
        Y = zeros(n, 1);
        X(1) = rand();
        Y(1) = rand();
        gamma_xy = 0.02;
        gamma_yx = 0.1;
        r_x = 3.7; r_y = 3.8;
        for i = 2:n
            X(i) = X(i - 1) * (r_x - r_x * X(i - 1) - gamma_xy * Y(i - 1));
            Y(i) = Y(i - 1) * (r_y - r_y * Y(i - 1) - gamma_yx * X(i - 1));
        end

    elseif mode == "RL"
        % n = 50
        alpha = 6;
        C = 2;
        dt = 0.01;
        X = zeros(n, 3);
        Y = zeros(n, 3);
        X(1, :) = (rand(1, 3) * 2 - 1) * 10;
        Y(1, :) = (rand(1, 3) * 2 - 1) * 10;
        
        cause_start = 1;
%         X(1,:) = ones(1,3) * 0.1;   Y(1,:) = ones(1,3) * 0.1;
        for i = 2:n
            X(i, 1) = X(i-1, 1) - alpha * (X(i-1, 2) + X(i-1, 3)) * dt;
            X(i, 2) = X(i-1, 2) + alpha * (X(i-1, 1) + 0.2 * X(i-1, 2)) * dt;
            X(i, 3) = X(i-1, 3) + alpha * (0.2 + X(i-1, 3) * (X(i-1, 1) - 5.7)) * dt;
            Y(i, 1) = Y(i-1, 1) + 10 * (-Y(i-1, 1) + Y(i-1, 2)) * dt;
            if i > cause_start
                Y(i, 2) = Y(i-1, 2) + (28 * Y(i-1, 1) - Y(i-1, 2) - Y(i-1, 1) * Y(i-1, 3) + C * X(i-cause_start, 2)^2) * dt;
            else
                Y(i, 2) = Y(i-1, 2) + (28 * Y(i-1, 1) - Y(i-1, 2) - Y(i-1, 1) * Y(i-1, 3)) * dt;
            end
            Y(i, 3) = Y(i-1, 3) + (Y(i-1, 1) * Y(i-1, 2) - 8/3 * Y(i-1, 3)) * dt;
            
            X(i,:) = X(i,:) + (rand(1,3)-0.5) * 0.1; Y(i,:) = Y(i,:) + (rand(1,3)-0.5) * 0.1;
        end
        
        X = X'; Y = Y';
    end

%     X = X(:, 101:end);
%     Y = Y(:, 101:end);
    X = X(101:end, :);
    Y = Y(101:end, :);

end