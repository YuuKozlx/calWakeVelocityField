function [] = splitF12(F12, each_group_size)
    % F12 is the data to be split
    % each_group_size is the size of each group. If it is not specified, it will be 1000 by default.
    if nargin < 2
        each_group_size = 1000;
    end

    group_numer = floor(length(F12) / group_size); % the number of groups

    for i = 1:group_numer

        if (i == group_numer) && (mod(length(F12), group_size) ~= 0)
            each_group_F12 = F12(i * each_group_size + 1:end, :);
            save(['F12', num2str(i * each_group_size + 1), '_', num2str(size(F12, 1)), '.mat'], 'each_group_F12');
        else
            each_group_F12 = F12((i - 1) * each_group_size + 1:i * each_group_size, :);
            save(['F12', num2str((i - 1) * each_group_size + 1), '_', num2str(i * each_group_size), '.mat'], 'each_group_F12');
        end

    end

end
