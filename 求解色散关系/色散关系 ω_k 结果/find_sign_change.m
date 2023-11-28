function result_indices = find_sign_change(A)

    % 初始化异号点的位置矩阵
    result_indices = [];

    % 找出前后两点异号的位置或等于0的位置
    for i = 1:size(A, 1)
        row = A(i, :);
        sign_changes = find((row(1:end - 1) .* row(2:end) < 0)); 
        zero_indices = find(abs(row) <= 1e-15);

        result = [sign_changes, zero_indices];
        result = sort(result);

        if length(result) >= 4 % 只取前4个
            result = result(end - 3:end);
        end

        result_indices = [result_indices; result];
    end
    result_indices = fliplr(result_indices);
end
