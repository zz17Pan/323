function [tx_array, rx_array] = init_arrays(params)
%INIT_ARRAYS 初始化发射和接收天线阵列
%   params: 系统参数结构体
%   tx_array: 发射阵列结构体
%   rx_array: 接收阵列结构体

% 发射阵列初始化
tx_array = struct();
tx_array.size = params.tx.array_size;
tx_array.pos = params.tx.pos;
tx_array.spacing = params.tx.spacing;

% 初始化发射阵列天线元素坐标
[tx_array.elements_pos, tx_array.num_elements] = create_array_elements(tx_array.size, tx_array.spacing, tx_array.pos);

% 接收阵列初始化
rx_array = struct();
rx_array.size = params.rx.array_size;
rx_array.pos = params.rx.init_pos;
rx_array.spacing = params.rx.spacing;

% 初始化接收阵列天线元素坐标
[rx_array.elements_pos, rx_array.num_elements] = create_array_elements(rx_array.size, rx_array.spacing, rx_array.pos);

% 检查接收阵列是否正确初始化
if rx_array.num_elements > 1
    first_pos = rx_array.elements_pos(1, :);
    last_pos = rx_array.elements_pos(end, :);
    
    if norm(first_pos - last_pos) < 1e-6
        fprintf('警告: 接收阵列元素未正确分布，正在重新初始化...\n');
        fprintf('阵列参数: 大小=%dx%d, 中心=[%.2f, %.2f, %.2f], 间距=%.2f\n', ...
            rx_array.size(1), rx_array.size(2), rx_array.pos(1), rx_array.pos(2), rx_array.pos(3), rx_array.spacing);
            
        % 强制重新初始化接收阵列元素
        [rx_array.elements_pos, ~] = create_array_elements(rx_array.size, rx_array.spacing, rx_array.pos);
        
        % 再次检查
        if norm(rx_array.elements_pos(1, :) - rx_array.elements_pos(end, :)) < 1e-6
            fprintf('错误: 接收阵列仍未正确初始化，尝试手动设置元素位置\n');
            
            % 手动设置阵元位置
            num_rows = rx_array.size(1);
            num_cols = rx_array.size(2);
            spacing = rx_array.spacing;
            center = rx_array.pos;
            
            % 计算阵列尺寸
            total_width = (num_cols - 1) * spacing;
            total_height = (num_rows - 1) * spacing;
            
            % 计算左上角位置
            left_top_x = center(1) - total_width/2;
            left_top_y = center(2) - total_height/2;
            left_top_z = center(3);
            
            % 手动填充阵元位置
            element_idx = 1;
            for row = 1:num_rows
                for col = 1:num_cols
                    x = left_top_x + (col-1) * spacing;
                    y = left_top_y + (row-1) * spacing;
                    z = left_top_z;
                    rx_array.elements_pos(element_idx, :) = [x, y, z];
                    element_idx = element_idx + 1;
                end
            end
        end
    end
end

% 打印接收阵列信息以便调试
fprintf('接收阵列初始位置: [%.2f, %.2f, %.2f]\n', rx_array.pos(1), rx_array.pos(2), rx_array.pos(3));
fprintf('接收阵列元素数量: %d\n', rx_array.num_elements);
fprintf('接收阵列第一个元素位置: [%.2f, %.2f, %.2f]\n', rx_array.elements_pos(1,1), rx_array.elements_pos(1,2), rx_array.elements_pos(1,3));
fprintf('接收阵列最后一个元素位置: [%.2f, %.2f, %.2f]\n', rx_array.elements_pos(end,1), rx_array.elements_pos(end,2), rx_array.elements_pos(end,3));

% 验证阵列元素是否正确分布
if rx_array.num_elements > 1
    if norm(rx_array.elements_pos(1, :) - rx_array.elements_pos(end, :)) < 1e-6
        fprintf('严重错误: 接收阵列元素仍未正确分布，请检查阵列初始化代码\n');
    else
        fprintf('接收阵列元素已正确分布，阵列尺寸: %.2f x %.2f m\n', ...
            norm(rx_array.elements_pos(1, 1:2) - rx_array.elements_pos(rx_array.size(2), 1:2)), ...
            norm(rx_array.elements_pos(1, 1:2) - rx_array.elements_pos((rx_array.size(1)-1)*rx_array.size(2)+1, 1:2)));
    end
end

end

function [elements_pos, num_elements] = create_array_elements(array_size, spacing, array_pos)
%CREATE_ARRAY_ELEMENTS 创建阵列中各天线元素的坐标
%   array_size: 阵列大小 [行, 列]
%   spacing: 阵元间距
%   array_pos: 阵列中心位置 [x, y, z]
%   elements_pos: 各天线元素的坐标 (Nx3矩阵, 每行为一个天线元素的[x,y,z]坐标)
%   num_elements: 天线元素总数

% 计算天线元素总数
num_rows = array_size(1);
num_cols = array_size(2);
num_elements = num_rows * num_cols;

% 初始化元素坐标矩阵
elements_pos = zeros(num_elements, 3);

% 计算阵列尺寸
total_width = (num_cols - 1) * spacing;
total_height = (num_rows - 1) * spacing;

% 阵列左上角坐标 (以阵列中心为参考点)
left_top_x = array_pos(1) - total_width/2;
left_top_y = array_pos(2) - total_height/2;
left_top_z = array_pos(3);

% 计算每个天线元素的绝对坐标
element_idx = 1;
for row = 1:num_rows
    for col = 1:num_cols
        % 计算当前元素坐标
        x = left_top_x + (col-1) * spacing;
        y = left_top_y + (row-1) * spacing;
        z = left_top_z;
        
        % 保存元素坐标
        elements_pos(element_idx, :) = [x, y, z];
        element_idx = element_idx + 1;
    end
end

end 