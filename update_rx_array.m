function rx_array = update_rx_array(rx_array, new_pos)
%UPDATE_RX_ARRAY 更新接收阵列位置
%   rx_array: 接收阵列结构体
%   new_pos: 新的阵列中心位置 [x, y, z]
%   返回更新后的接收阵列结构体

% 检查接收阵列元素是否正确展开
% 如果所有元素都在同一个位置，则重新初始化阵列
if rx_array.num_elements > 1
    first_pos = rx_array.elements_pos(1, :);
    last_pos = rx_array.elements_pos(end, :);
    
    if norm(first_pos - last_pos) < 1e-6
        fprintf('警告: 检测到接收阵列元素未正确展开，调用create_array_elements强制重新初始化\n');
        
        % 当前阵列位置使用old_pos，因为尚未更新
        old_pos = rx_array.pos;
        
        % 调用create_array_elements重新初始化阵元
        [rx_array.elements_pos, ~] = create_array_elements(rx_array.size, rx_array.spacing, old_pos);
        
        % 验证阵元是否已正确初始化
        if norm(rx_array.elements_pos(1, :) - rx_array.elements_pos(end, :)) > 1e-6
            fprintf('  阵列已重新初始化 - 第一个元素: [%.2f, %.2f, %.2f], 最后一个元素: [%.2f, %.2f, %.2f]\n', ...
                rx_array.elements_pos(1,1), rx_array.elements_pos(1,2), rx_array.elements_pos(1,3), ...
                rx_array.elements_pos(end,1), rx_array.elements_pos(end,2), rx_array.elements_pos(end,3));
        else
            fprintf('错误: 接收阵列重新初始化失败，手动设置阵列元素位置\n');
            
            % 手动创建阵列元素位置
            [rx_array.elements_pos, ~] = manual_create_array_elements(rx_array.size, rx_array.spacing, old_pos);
        end
    end
end

% 计算位置偏移量
pos_offset = new_pos - rx_array.pos;

% 更新阵列中心位置
rx_array.pos = new_pos;

% 更新每个天线元素的位置
for i = 1:rx_array.num_elements
    rx_array.elements_pos(i, :) = rx_array.elements_pos(i, :) + pos_offset;
end

% 打印更新后的阵列信息
fprintf('更新后的接收阵列位置: [%.2f, %.2f, %.2f]\n', rx_array.pos(1), rx_array.pos(2), rx_array.pos(3));
fprintf('第一个元素位置: [%.2f, %.2f, %.2f], 最后一个元素位置: [%.2f, %.2f, %.2f]\n', ...
    rx_array.elements_pos(1,1), rx_array.elements_pos(1,2), rx_array.elements_pos(1,3), ...
    rx_array.elements_pos(end,1), rx_array.elements_pos(end,2), rx_array.elements_pos(end,3));

% 确认元素位置已正确更新
if rx_array.num_elements > 1 && norm(rx_array.elements_pos(1, :) - rx_array.elements_pos(end, :)) < 1e-6
    fprintf('严重错误: 更新后接收阵列元素仍然集中在同一点!\n');
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

function [elements_pos, num_elements] = manual_create_array_elements(array_size, spacing, array_pos)
%MANUAL_CREATE_ARRAY_ELEMENTS 手动创建阵列中各天线元素的坐标
%   这是create_array_elements的备份函数，用于当正常方法失效时

% 计算天线元素总数
num_rows = array_size(1);
num_cols = array_size(2);
num_elements = num_rows * num_cols;

% 初始化元素坐标矩阵
elements_pos = zeros(num_elements, 3);

% 强制设置阵列间距
if spacing < 0.01
    fprintf('  警告: 阵元间距过小 (%.4f)，强制设置为0.1米\n', spacing);
    spacing = 0.1;
end

% 计算阵列在每个维度的总宽度
width_x = (num_cols - 1) * spacing;
width_y = (num_rows - 1) * spacing;

% 明确计算左上角起始位置
start_x = array_pos(1) - width_x/2;
start_y = array_pos(2) - width_y/2;
start_z = array_pos(3);

fprintf('  手动创建阵列: %dx%d, 中心=[%.2f,%.2f,%.2f], 间距=%.2f\n', ...
    num_rows, num_cols, array_pos(1), array_pos(2), array_pos(3), spacing);
fprintf('  阵列左上角: [%.2f,%.2f,%.2f], 尺寸: %.2f x %.2f\n', ...
    start_x, start_y, start_z, width_x, width_y);

% 填充阵元位置
element_idx = 1;
for row = 1:num_rows
    for col = 1:num_cols
        % 明确计算每个元素的绝对位置
        x = start_x + (col-1) * spacing;
        y = start_y + (row-1) * spacing;
        z = start_z;
        
        % 保存到结果矩阵
        elements_pos(element_idx, :) = [x, y, z];
        
        % 打印调试信息
        if element_idx == 1 || element_idx == num_elements
            fprintf('  元素 %d 位置: [%.2f, %.2f, %.2f]\n', element_idx, x, y, z);
        end
        
        element_idx = element_idx + 1;
    end
end

end 