function X_3d = convert_to_3d_array(X)
% CONVERT_TO_3D_ARRAY Converts input to a 3D array
%
% Input:
%   X - Either a cell array of 2D arrays, or a 2D array
%
% Output:
%   X_3d - A 3D array of size [N_arrays, rows, cols]
%
% Usage:
%   X_3d = convert_to_3d_array(X);
%
% Examples:
%   % For a cell array of 2D arrays
%   X = {rand(5,10), rand(5,10), rand(5,10)};
%   X_3d = convert_to_3d_array(X); % Returns a 3x5x10 array
%
%   % For a 2D array
%   Y = rand(5,10);
%   Y_3d = convert_to_3d_array(Y); % Returns a 1x5x10 array

% Check if X is a cell array
if iscell(X)
    % Check if all elements are 2D arrays
    all_2d = true;
    for i = 1:length(X)
        if ~isnumeric(X{i}) || ndims(X{i}) ~= 2
            all_2d = false;
            break;
        end
    end

    if all_2d
        % Check if all arrays have the same dimensions
        % dims = cellfun(@size, X, 'UniformOutput', false);
        % if ~all(cellfun(@(d) isequal(d, dims{1}), dims))
        %     error('All 2D arrays in the cell array must have the same dimensions');
        % end

        % Get dimensions
        [rows, cols] = size(X{1});
        num_arrays = length(X);

        % Create 3D array: N_arrays x rows x cols
        X_3d = zeros(num_arrays, rows, cols);

        % Fill the 3D array
        for i = 1:num_arrays
            X_3d(i,:,:) = X{i};
        end
    else
        error('Not all elements in the cell array are 2D numeric arrays');
    end
elseif isnumeric(X) && ismatrix(X)
    % If X is a 2D array, treat it as a single array: 1 x rows x cols
    X_3d = permute(X, [3 1 2]);
elseif isnumeric(X) && ndims(X) == 3
    % If X is a 3D array, assume it is already N_arrays x rows x cols
    X_3d = X;
else
    error('Input must be either a cell array of 2D numeric arrays or a 2D/3D numeric array');
end

end
