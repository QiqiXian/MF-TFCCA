% compute the year-to-year growth rates of a time series data, p data points per year

function [yoy] = YOY(data, p, year_start, year_end)
    % data: N * 3 (year, month, value)

    % find the starting and ending indices
    n = size(data, 1);
    start_idx = 0;
    end_idx = 0;
    for i = 1:n
        if data(i, 1) == year_start
            start_idx = i;
            break;
        end
    end
    for i = n:-1:1
        if data(i, 1) == year_end
            end_idx = i;
            break;
        end
    end

    data = data(start_idx:end_idx, :);
    data = reshape(data(:, 3), p, [])';
%     data = reshape(data, p, [])';
    nyears = year_end - year_start + 1;
%     nyears = size(data, 1);

    % year-to-year growth rates
    yoy = zeros(nyears - 1, p);
    for i = 1:nyears - 1
        for j = 1:p
            yoy(i, j) = ((data(i + 1, j)) / data(i, j) - 1) * 100;
        end
    end

end