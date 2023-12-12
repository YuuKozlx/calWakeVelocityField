function result_indices = find_sign_change(A, zeronum)
    % 输入参数: 矩阵A zeronum:零点的个数
    % 输出参数: 矩阵A中每一行前后两个异号的位置或等于0的位置
    % 这个函数的作用:
    % 根据王宏伟博士论文中 公式(4-61) 「 F12 = 0 」求解零点的位置
    % 矩阵A中每一行前后两个异号的位置或等于0的位置
    % 但是输出零点的位置是从后往前计算的
    % 例如: [1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 0, 0, 0]
    % 输出: [10 11 12 13]

    if nargin < 2
        zeronum = 3;
    end

    % 初始化异号点的位置矩阵
    result_indices = zeros(size(A, 1), zeronum);

    % 找出前后两点异号的位置或等于0的位置
    for i = 1:size(A, 1)
        row = A(i, :);
        sign_changes = find((row(1:end - 1) .* row(2:end) < 0)); % 前后两点异号的位置
        zero_indices = find(abs(row) <= 1e-15); % 等于0的位置,由于计算误差,可选取绝对值小于1e-15的数,认为是0

        result = [sign_changes, zero_indices];
        result = sort(result);

        if length(result) >= zeronum % 只取前zeronum个
            result = result(end - zeronum+1:end);
        end

        result_indices(i, :) = result;
    end

    result_indices = fliplr(result_indices);
end