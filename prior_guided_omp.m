function [estimated_range, estimated_azimuth, estimated_elevation] = prior_guided_omp(rx_signal, tx_array, rx_array, prior, prior_covariance, params)
%PRIOR_GUIDED_OMP 使用先验信息引导的OMP算法进行参数估计
%   使用正交匹配追踪算法进行稀疏重建，结合先验信息提高估计精度
%   优化版本：减少计算复杂度，提高性能和速度

% 设置网格参数 - 减小网格密度以降低计算复杂度
MAX_GRID_POINTS = 7;  % 最大网格点数（由10降为7）
MIN_GRID_POINTS = 3;  % 最小网格点数（由5降为3）

% 提取先验信息
prior_range = prior.range;
prior_azimuth = prior.azimuth;
prior_elevation = prior.elevation;

% 提取协方差信息
range_var = prior_covariance(1,1);
azimuth_var = prior_covariance(2,2);
elevation_var = prior_covariance(3,3);

% 计算标准差
range_std = sqrt(range_var);
azimuth_std = sqrt(azimuth_var);
elevation_std = sqrt(elevation_var);

% 计算搜索范围 - 更严格的范围以减少计算量
% 根据先验信息和标准差设置合理的搜索范围
% 确保最小距离不会太小，避免出现bin超出范围的问题
min_range = max(10, prior_range - 2.0 * range_std);  % 缩小为2.0倍标准差（原为3.0），且不小于10米
max_range = min(3000, prior_range + 2.0 * range_std); % 缩小为2.0倍标准差（原为3.0）

% 确保最小和最大距离差别至少有10米，避免搜索范围过小
if max_range - min_range < 10
    center_range = (max_range + min_range) / 2;
    min_range = center_range - 5;
    max_range = center_range + 5;
end

% 角度搜索范围更加严格
min_azimuth = max(-180, prior_azimuth - 1.5 * azimuth_std);  % 缩小为1.5倍标准差（原为2.0）
max_azimuth = min(180, prior_azimuth + 1.5 * azimuth_std);   % 缩小为1.5倍标准差（原为2.0）
min_elevation = max(-90, prior_elevation - 1.5 * elevation_std);  % 缩小为1.5倍标准差（原为2.0）
max_elevation = min(90, prior_elevation + 1.5 * elevation_std);   % 缩小为1.5倍标准差（原为2.0）

% 计算理论上的最大距离，避免过大的搜索范围
c = 3e8;  % 光速
max_theoretical_range = params.c * params.fmcw.fs / (2 * params.fmcw.mu * 2);
if max_range > max_theoretical_range * 0.95
    fprintf('警告: 搜索距离范围 (%.2f m) 接近理论最大距离 (%.2f m), 调整为理论最大值的90%%\n', ...
            max_range, max_theoretical_range);
    max_range = max_theoretical_range * 0.9;
end

% 适应性确定网格点数，减少计算量
range_span = max_range - min_range;
azimuth_span = max_azimuth - min_azimuth;
elevation_span = max_elevation - min_elevation;

% 根据搜索范围自适应调整网格点数
range_grid_points = min(MAX_GRID_POINTS, max(MIN_GRID_POINTS, ceil(range_span / 10)));  
azimuth_grid_points = min(MAX_GRID_POINTS, max(MIN_GRID_POINTS, ceil(azimuth_span / 5)));
elevation_grid_points = min(MAX_GRID_POINTS, max(MIN_GRID_POINTS, ceil(elevation_span / 5)));

% 创建搜索网格
range_grid = linspace(min_range, max_range, range_grid_points);
azimuth_grid = linspace(min_azimuth, max_azimuth, azimuth_grid_points);
elevation_grid = linspace(min_elevation, max_elevation, elevation_grid_points);

fprintf('OMP优化搜索范围: 距离=[%.2f, %.2f] m (%d点), 方位角=[%.2f, %.2f]° (%d点), 俯仰角=[%.2f, %.2f]° (%d点)\n', ...
        min_range, max_range, range_grid_points, ...
        min_azimuth, max_azimuth, azimuth_grid_points, ...
        min_elevation, max_elevation, elevation_grid_points);

% 提取信号参数
[Ns, num_chirps, num_rx] = size(rx_signal);
num_tx = tx_array.num_elements;

% 性能优化：只处理与最大距离对应的FFT点数，而不是全部采样点
% 计算需要的FFT点数
% 计算理论上的最小有效距离对应的bin
min_valid_range = 1; % 至少1米以上有效
min_valid_bin = ceil(min_valid_range * 2 * params.fmcw.mu / params.c * params.fmcw.fs);

% 确保最小距离bin不会导致警告
if min_range < min_valid_range
    min_range = min_valid_range;
    fprintf('已自动调整最小搜索距离为 %d 米，避免bin超出范围\n', min_valid_range);
end

max_fft_bin = ceil(max_range * 2 * params.fmcw.mu / params.c * params.fmcw.fs);
max_fft_bin = min(max_fft_bin, Ns);  % 确保不超过采样点数
min_fft_bin = ceil(min_range * 2 * params.fmcw.mu / params.c * params.fmcw.fs);
min_fft_bin = max(1, min_fft_bin);  % 确保至少从第1个bin开始

% 检查最小和最大bin的有效性
if min_fft_bin >= Ns
    fprintf('警告: 最小距离bin(%d)超出采样点数(%d)，调整为1\n', min_fft_bin, Ns);
    min_fft_bin = 1;
end

if max_fft_bin <= min_fft_bin
    fprintf('最大距离bin(%d)小于等于最小距离bin(%d)，自动调整范围\n', max_fft_bin, min_fft_bin);
    max_fft_bin = min(min_fft_bin + 100, Ns);
end

% 为了效率，只计算感兴趣范围内的FFT
% 选择更安全的FFT长度
fft_length = 2^nextpow2(Ns);  % 使用基于采样点数的幂次，而非基于bin索引
if fft_length > 2^15  % 限制最大长度为32768
    fft_length = 2^15;
end

% 确保索引不超过FFT长度
if max_fft_bin >= fft_length
    fprintf('警告: 最大FFT bin索引(%d)超过FFT长度(%d)，进行调整\n', max_fft_bin, fft_length);
    max_fft_bin = fft_length - 1;
end

range_bin_indices = round(linspace(min_fft_bin, max_fft_bin, range_grid_points));
% 再次检查所有索引是否在有效范围内
range_bin_indices = max(1, min(range_bin_indices, fft_length-1));

% 性能优化：使用较少的chirps来降低计算量
max_chirps_to_use = min(num_chirps, 16);  % 最多使用16个chirps
chirp_indices = round(linspace(1, num_chirps, max_chirps_to_use));
num_chirps_used = length(chirp_indices);

fprintf('性能优化: 使用 %d/%d chirps, FFT长度=%d, 距离bin范围=[%d, %d]\n', ...
        num_chirps_used, num_chirps, fft_length, min_fft_bin, max_fft_bin);

% 创建局部距离FFT和全体信号模型
local_fft = zeros(length(range_bin_indices), num_chirps_used, num_rx);

% 只对感兴趣的距离范围执行FFT，而不是全部
for chirp_idx = 1:num_chirps_used
    actual_chirp = chirp_indices(chirp_idx);
    for rx_idx = 1:num_rx
        % 对当前chirp和接收通道执行FFT
        % 安全处理：确保FFT长度不大于信号长度
        safe_fft_length = min(fft_length, size(rx_signal, 1));
        chirp_fft = fft(rx_signal(:, actual_chirp, rx_idx), safe_fft_length);
        
        % 再次验证索引范围
        valid_indices = find(range_bin_indices > 0 & range_bin_indices <= length(chirp_fft));
        
        % 只保留有效的距离bin
        for i = 1:length(valid_indices)
            idx = valid_indices(i);
            bin_idx = range_bin_indices(idx);
            if bin_idx <= length(chirp_fft)
                local_fft(idx, chirp_idx, rx_idx) = chirp_fft(bin_idx);
            end
        end
    end
end

% 初始化OMP算法
max_iter = params.omp.max_iter;
residual = local_fft;  % 初始残差
estimated_signal = zeros(size(local_fft));  % 估计信号
selected_indices = [];  % 已选择的字典原子索引

% 预计算先验权重 - 提高性能
% 根据先验信息和协方差计算权重
range_prior_weights = exp(-(range_grid - prior_range).^2 / (2 * range_var));
azimuth_prior_weights = exp(-(azimuth_grid - prior_azimuth).^2 / (2 * azimuth_var));
elevation_prior_weights = exp(-(elevation_grid - prior_elevation).^2 / (2 * elevation_var));

% 归一化权重
range_prior_weights = range_prior_weights / sum(range_prior_weights);
azimuth_prior_weights = azimuth_prior_weights / sum(azimuth_prior_weights);
elevation_prior_weights = elevation_prior_weights / sum(elevation_prior_weights);

% 使用更高效的OMP算法 - 混合全局和局部搜索
for iter = 1:max_iter
    fprintf('OMP迭代 %d/%d...\n', iter, max_iter);
    
    % 变量初始化
    max_correlation = -inf;
    best_r_idx = 1; % 初始化为有效索引1，而不是0
    best_az_idx = 1; % 初始化为有效索引1，而不是0
    best_el_idx = 1; % 初始化为有效索引1，而不是0
    
    if iter == 1
        % 第一次迭代: 全局网格搜索
        for r_idx = 1:length(range_grid)
            range_val = range_grid(r_idx);
            range_weight = range_prior_weights(r_idx);
            
            for az_idx = 1:length(azimuth_grid)
                azimuth_val = azimuth_grid(az_idx);
                azimuth_weight = azimuth_prior_weights(az_idx);
                
                for el_idx = 1:length(elevation_grid)
                    elevation_val = elevation_grid(el_idx);
                    elevation_weight = elevation_prior_weights(el_idx);
                    
                    % 计算当前参数的字典原子
                    atom = generate_steering_vector(range_val, azimuth_val, elevation_val, ...
                                                  range_bin_indices, chirp_indices, tx_array, rx_array, params);
                    
                    % 计算与残差的相关性，加入先验权重
                    correlation = abs(sum(residual .* conj(atom), 'all'));
                    weighted_correlation = correlation * range_weight * azimuth_weight * elevation_weight;
                    
                    % 更新最大相关性
                    if weighted_correlation > max_correlation
                        max_correlation = weighted_correlation;
                        best_r_idx = r_idx;
                        best_az_idx = az_idx;
                        best_el_idx = el_idx;
                    end
                end
            end
        end
    else
        % 后续迭代：局部搜索 - 在上一次最优点附近进行更细致的搜索
        % 确保索引有效
        best_r_idx = max(1, min(best_r_idx, length(range_grid)));
        best_az_idx = max(1, min(best_az_idx, length(azimuth_grid)));
        best_el_idx = max(1, min(best_el_idx, length(elevation_grid)));
        
        % 生成局部精细网格 - 仅在最优点附近搜索
        local_r_range = linspace(max(min_range, range_grid(best_r_idx) - range_span/range_grid_points), ...
                                 min(max_range, range_grid(best_r_idx) + range_span/range_grid_points), ...
                                 5);
        local_az_range = linspace(max(min_azimuth, azimuth_grid(best_az_idx) - azimuth_span/azimuth_grid_points), ...
                                  min(max_azimuth, azimuth_grid(best_az_idx) + azimuth_span/azimuth_grid_points), ...
                                  5);
        local_el_range = linspace(max(min_elevation, elevation_grid(best_el_idx) - elevation_span/elevation_grid_points), ...
                                  min(max_elevation, elevation_grid(best_el_idx) + elevation_span/elevation_grid_points), ...
                                  5);
                                  
        % 在局部网格上搜索最优点
        % 初始化局部搜索的最优值
        best_r_val = range_grid(best_r_idx);
        best_az_val = azimuth_grid(best_az_idx);
        best_el_val = elevation_grid(best_el_idx);
        
        for r_val = local_r_range
            % 计算当前距离的先验权重
            r_weight = exp(-(r_val - prior_range)^2 / (2 * range_var));
            
            for az_val = local_az_range
                % 计算当前方位角的先验权重
                az_weight = exp(-(az_val - prior_azimuth)^2 / (2 * azimuth_var));
                
                for el_val = local_el_range
                    % 计算当前俯仰角的先验权重
                    el_weight = exp(-(el_val - prior_elevation)^2 / (2 * elevation_var));
                    
                    % 生成字典原子
                    atom = generate_steering_vector(r_val, az_val, el_val, ...
                                                   range_bin_indices, chirp_indices, tx_array, rx_array, params);
                    
                    % 计算相关性并带入先验权重
                    correlation = abs(sum(residual .* conj(atom), 'all'));
                    weighted_correlation = correlation * r_weight * az_weight * el_weight;
                    
                    % 更新最优值
                    if weighted_correlation > max_correlation
                        max_correlation = weighted_correlation;
                        best_r_val = r_val;
                        best_az_val = az_val;
                        best_el_val = el_val;
                    end
                end
            end
        end
        
        % 检查是否找到了更好的局部最优解
        if max_correlation > -inf
            % 更新最优值到网格中对应的索引
            [~, best_r_idx] = min(abs(range_grid - best_r_val));
            [~, best_az_idx] = min(abs(azimuth_grid - best_az_val));
            [~, best_el_idx] = min(abs(elevation_grid - best_el_val));
        else
            % 如果没有找到更好的局部最优解，保持原来的索引不变
            fprintf('  警告: 局部搜索未找到更好的解，使用上一次迭代的结果\n');
        end
    end
    
    % 获取最优参数
    best_range = range_grid(best_r_idx);
    best_azimuth = azimuth_grid(best_az_idx);
    best_elevation = elevation_grid(best_el_idx);
    
    fprintf('  找到最优参数: 距离=%.2f m, 方位角=%.2f°, 俯仰角=%.2f°\n', ...
            best_range, best_azimuth, best_elevation);
    
    % 生成最优字典原子
    best_atom = generate_steering_vector(best_range, best_azimuth, best_elevation, ...
                                       range_bin_indices, chirp_indices, tx_array, rx_array, params);
    
    % 归一化
    best_atom = best_atom / norm(best_atom(:));
    
    % 计算投影系数
    projection_coef = sum(residual .* conj(best_atom), 'all');
    
    % 更新估计信号
    estimated_signal = estimated_signal + projection_coef * best_atom;
    
    % 更新残差
    residual = residual - projection_coef * best_atom;
    
    % 记录选择的参数
    selected_indices = [selected_indices; best_r_idx, best_az_idx, best_el_idx];
    
    % 检查残差能量
    residual_energy = sum(abs(residual).^2, 'all');
    fprintf('  剩余残差能量: %.2e\n', residual_energy);
    
    % 残差能量足够小时终止迭代
    if residual_energy < 1e-4 * sum(abs(local_fft).^2, 'all')
        fprintf('  残差能量足够小，提前终止迭代\n');
        break;
    end
end

% 最终估计 - 根据迭代次数调整权重，最后一次迭代的权重更高
num_selected = size(selected_indices, 1);
if num_selected > 0
    % 给后期迭代的结果更高的权重
    iteration_weights = (1:num_selected) / sum(1:num_selected);
    
    % 初始化加权平均结果
    weighted_range = 0;
    weighted_azimuth = 0;
    weighted_elevation = 0;
    
    for i = 1:num_selected
        r_idx = selected_indices(i, 1);
        az_idx = selected_indices(i, 2);
        el_idx = selected_indices(i, 3);
        
        weighted_range = weighted_range + range_grid(r_idx) * iteration_weights(i);
        weighted_azimuth = weighted_azimuth + azimuth_grid(az_idx) * iteration_weights(i);
        weighted_elevation = weighted_elevation + elevation_grid(el_idx) * iteration_weights(i);
    end
    
    % 最终估计结果
    estimated_range = weighted_range;
    estimated_azimuth = weighted_azimuth;
    estimated_elevation = weighted_elevation;
    
    % 合理性检查 - 确保估计结果在物理上有意义
    % 距离估计合理性检查
    if estimated_range > max_theoretical_range
        warning('估计距离 (%.2f m) 超过理论最大距离 (%.2f m)，调整为理论最大值', estimated_range, max_theoretical_range);
        estimated_range = max_theoretical_range * 0.95;
    end
    
    if estimated_range < 1
        warning('估计距离 (%.2f m) 太小，调整为最小值 1米', estimated_range);
        estimated_range = 1;
    end
    
    % 角度估计合理性检查
    estimated_azimuth = mod(estimated_azimuth + 180, 360) - 180;  % 规范化到[-180,180]
    estimated_elevation = max(-90, min(90, estimated_elevation));  % 规范化到[-90,90]
    
    % 如果与先验差异过大，进行警告并适当调整
    if abs(estimated_range - prior_range) > 3 * range_std
        warning('估计距离 (%.2f m) 与先验 (%.2f m) 差异过大，调整为先验与估计的加权平均', estimated_range, prior_range);
        estimated_range = 0.7 * estimated_range + 0.3 * prior_range;
    end
    
    if abs(estimated_azimuth - prior_azimuth) > 3 * azimuth_std
        warning('估计方位角 (%.2f°) 与先验 (%.2f°) 差异过大，调整为先验与估计的加权平均', estimated_azimuth, prior_azimuth);
        estimated_azimuth = 0.7 * estimated_azimuth + 0.3 * prior_azimuth;
    end
    
    if abs(estimated_elevation - prior_elevation) > 3 * elevation_std
        warning('估计俯仰角 (%.2f°) 与先验 (%.2f°) 差异过大，调整为先验与估计的加权平均', estimated_elevation, prior_elevation);
        estimated_elevation = 0.7 * estimated_elevation + 0.3 * prior_elevation;
    end
else
    % 如果没有选择任何字典原子，返回先验值
    warning('OMP算法未能选择有效的字典原子，使用先验值作为估计结果');
    estimated_range = prior_range;
    estimated_azimuth = prior_azimuth;
    estimated_elevation = prior_elevation;
end

fprintf('OMP最终估计: 距离=%.2f m, 方位角=%.2f°, 俯仰角=%.2f°\n', ...
        estimated_range, estimated_azimuth, estimated_elevation);

end

% 辅助函数: 生成特定参数的字典原子
function steering_vector = generate_steering_vector(range, azimuth, elevation, range_bins, chirp_indices, tx_array, rx_array, params)
    % 生成指定参数下的导向矢量
    % range: 距离(米)
    % azimuth: 方位角(度)
    % elevation: 俯仰角(度)
    % range_bins: 距离bin索引向量
    % chirp_indices: chirp索引向量
    % tx_array, rx_array: 发射和接收阵列结构体
    % params: 系统参数
    
    % 波长
    lambda = params.c / params.fc;
    
    % 获取收发阵列的位置信息
    if isfield(tx_array, 'elements_pos')
        tx_positions = tx_array.elements_pos;
    else
        error('发射阵列缺少元素位置信息(elements_pos)');
    end
    
    if isfield(rx_array, 'elements_pos')
        rx_positions = rx_array.elements_pos;
    else
        error('接收阵列缺少元素位置信息(elements_pos)');
    end
    
    % 获取接收阵列阵元数
    num_rx = size(rx_positions, 1);
    
    % 计算方向向量
    sin_az = sind(azimuth);
    cos_az = cosd(azimuth);
    sin_el = sind(elevation);
    cos_el = cosd(elevation);
    
    direction = [cos_az * cos_el; sin_az * cos_el; sin_el];
    
    % 初始化透射矢量
    num_range_bins = length(range_bins);
    num_chirps = length(chirp_indices);
    steering_vector = zeros(num_range_bins, num_chirps, num_rx);
    
    % 计算基于距离的相位延迟
    range_phase = exp(-1j * 4 * pi * range * params.fmcw.mu / params.c * (1:num_range_bins) / params.fmcw.fs);
    
    % 为每个接收器、发射器和chirp计算相位
    for rx_idx = 1:num_rx
        % elements_pos是Nx3矩阵，每行代表一个元素的[x,y,z]坐标
        rx_pos = rx_positions(rx_idx, :)';  % 转置为列向量便于dot运算
        
        for tx_idx = 1:1  % 简化：只考虑第一个发射器
            tx_pos = tx_positions(tx_idx, :)';  % 转置为列向量
            
            % 计算到当前发射器和接收器的总路径长度
            total_path = range;  % 简化：使用提供的距离而不是计算精确路径
            
            % 计算角度相关的相位延迟 - 确保所有向量都是列向量
            phase_delay = 2 * pi / lambda * (dot(rx_pos, direction) + dot(tx_pos, direction));
            
            % 为每个chirp和距离bin计算相位
            for chirp_idx = 1:num_chirps
                actual_chirp = chirp_indices(chirp_idx);
                % 修正这里：添加确保参数存在的检查
                if isfield(params.fmcw, 'chirps_per_frame')
                    doppler_phase = exp(-1j * 2 * pi * 2 * range / lambda * actual_chirp / params.fmcw.chirps_per_frame);
                else
                    % 使用num_chirps作为替代方案
                    doppler_phase = exp(-1j * 2 * pi * 2 * range / lambda * actual_chirp / params.fmcw.num_chirps);
                end
                
                for range_idx = 1:num_range_bins
                    bin_phase = range_phase(range_idx);
                    steering_vector(range_idx, chirp_idx, rx_idx) = ...
                        steering_vector(range_idx, chirp_idx, rx_idx) + bin_phase * doppler_phase * exp(1j * phase_delay);
                end
            end
        end
    end
    
    % 归一化
    steering_vector = steering_vector / sqrt(sum(abs(steering_vector).^2, 'all'));
end 