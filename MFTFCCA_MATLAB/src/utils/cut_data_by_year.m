function [data] = cut_data_by_year(data, year_start, year_end)
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
end
