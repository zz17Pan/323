function [estimated_azimuth, estimated_elevation] = music_angle_estimation(rx_signal, range_idx, doppler_idx, tx_array, rx_array, params, prior_azimuth, prior_elevation, az_std, el_std)
%MUSIC_ANGLE_ESTIMATION 使用MUSIC算法进行角度估计
%   优化的MUSIC实现，减少搜索空间，提高计算性能
%   rx_signal: 接收信号，尺寸为 [samples, chirps, rx_elements]
%   range_idx: 目标对应的距离索引
%   doppler_idx: 目标对应的多普勒索引（如果为空，使用零多普勒）
%   tx_array, rx_array: 发射和接收阵列结构体
%   params: 系统参数结构体
%   prior_azimuth, prior_elevation: 先验角度估计（可选）
%   az_std, el_std: 先验角度的标准差（可选）

% 处理可选参数
if nargin < 7 || isempty(prior_azimuth)
    prior_azimuth = 0;  % 默认值
end

if nargin < 8 || isempty(prior_elevation)
    prior_elevation = 0;  % 默认值
end

if nargin < 9 || isempty(az_std)
    az_std = 45;  % 默认搜索范围较大
end

if nargin < 10 || isempty(el_std)
    el_std = 45;  % 默认搜索范围较大
end

% 性能优化：只提取我们需要的范围bin数据
range_bin_signal = squeeze(rx_signal(range_idx, :, :));  % [chirps, rx_elements]

% 性能优化：限制处理的chirp数量以降低计算量
max_chirps = 32;  % 最多使用32个chirps
if size(range_bin_signal, 1) > max_chirps
    % 均匀采样chirps
    chirp_indices = round(linspace(1, size(range_bin_signal, 1), max_chirps));
    range_bin_signal = range_bin_signal(chirp_indices, :);
    fprintf('MUSIC性能优化: 从%d个chirps降至%d个\n', size(rx_signal, 2), max_chirps);
end

% 计算协方差矩阵并添加对角加载以提高数值稳定性
R = range_bin_signal' * range_bin_signal / size(range_bin_signal, 1);
R = R + eye(size(R)) * (1e-6 * trace(R) / size(R, 1));  % 对角加载

% 处理先验信息
has_prior = true;  % 现在我们总是有先验信息，即使是默认值
% 根据先验信息设置搜索范围
if abs(prior_azimuth) < 1e-6 && abs(prior_elevation) < 1e-6
    % 如果先验值接近零（默认值），使用全范围搜索
    az_min = params.music.az_range(1);
    az_max = params.music.az_range(2);
    el_min = params.music.el_range(1);
    el_max = params.music.el_range(2);
else
    % 根据先验信息和标准差限制搜索范围
    az_min = max(params.music.az_range(1), prior_azimuth - 3*az_std);
    az_max = min(params.music.az_range(2), prior_azimuth + 3*az_std);
    el_min = max(params.music.el_range(1), prior_elevation - 3*el_std);
    el_max = min(params.music.el_range(2), prior_elevation + 3*el_std);
end

% 计算特征分解
[V, D] = eig(R);
[eigvals, idx] = sort(diag(D), 'descend');
V = V(:, idx);

% 确定信号子空间维度
% 自适应确定信号子空间维度
eigvals_norm = eigvals / max(eigvals);
signal_idx = find(eigvals_norm > 0.1, 1, 'last');  % 使用10%阈值
if isempty(signal_idx) || signal_idx < 1
    signal_idx = 1;  % 至少有一个信号
end
noise_subspace = V(:, signal_idx+1:end);

% 两步搜索策略：先粗后细
% 1. 粗搜索
coarse_az_steps = 15;  % 粗搜索方位角步长
coarse_el_steps = 10;  % 粗搜索俯仰角步长
coarse_az_grid = linspace(az_min, az_max, coarse_az_steps);
coarse_el_grid = linspace(el_min, el_max, coarse_el_steps);

% 初始化粗搜索结果
coarse_spectrum = zeros(coarse_az_steps, coarse_el_steps);

% 计算粗搜索MUSIC谱
for az_idx = 1:coarse_az_steps
    az = coarse_az_grid(az_idx);
    for el_idx = 1:coarse_el_steps
        el = coarse_el_grid(el_idx);
        a = steering_vector(az, el, rx_array.elements_pos, params);
        coarse_spectrum(az_idx, el_idx) = 1 / (a' * (noise_subspace * noise_subspace') * a);
    end
end

% 找到粗搜索最大值位置
[~, max_idx] = max(coarse_spectrum(:));
[max_az_idx, max_el_idx] = ind2sub(size(coarse_spectrum), max_idx);
coarse_az = coarse_az_grid(max_az_idx);
coarse_el = coarse_el_grid(max_el_idx);

% 2. 精细搜索 - 围绕粗搜索最大值进行
fine_range = 1.5 * max(coarse_az_steps, coarse_el_steps) / coarse_az_steps;  % 精细搜索范围
fine_az_range = [max(az_min, coarse_az - fine_range), min(az_max, coarse_az + fine_range)];
fine_el_range = [max(el_min, coarse_el - fine_range), min(el_max, coarse_el + fine_range)];
fine_az_steps = 30;  % 精细搜索步数
fine_el_steps = 20;  % 精细搜索步数
fine_az_grid = linspace(fine_az_range(1), fine_az_range(2), fine_az_steps);
fine_el_grid = linspace(fine_el_range(1), fine_el_range(2), fine_el_steps);

% 初始化精细搜索谱
fine_spectrum = zeros(fine_az_steps, fine_el_steps);

% 性能优化：预计算用于向量化的矩阵
az_meshgrid = repmat(fine_az_grid(:), 1, fine_el_steps);
el_meshgrid = repmat(fine_el_grid, fine_az_steps, 1);
all_steering_vectors = zeros(size(rx_array.elements_pos, 1), fine_az_steps * fine_el_steps);

% 批量计算所有方向的导向矢量
for i = 1:(fine_az_steps * fine_el_steps)
    az = az_meshgrid(i);
    el = el_meshgrid(i);
    all_steering_vectors(:, i) = steering_vector(az, el, rx_array.elements_pos, params);
end

% 批量计算MUSIC谱
P = noise_subspace * noise_subspace';
spectrum_values = zeros(fine_az_steps * fine_el_steps, 1);

for i = 1:(fine_az_steps * fine_el_steps)
    a = all_steering_vectors(:, i);
    spectrum_values(i) = 1 / real(a' * P * a);  % 取实部以避免数值误差
end

% 重塑谱到网格形式
fine_spectrum = reshape(spectrum_values, fine_az_steps, fine_el_steps);

% 找出精细谱的最大值
[peak_value, max_idx] = max(fine_spectrum(:));
[max_az_idx, max_el_idx] = ind2sub(size(fine_spectrum), max_idx);
estimated_azimuth = fine_az_grid(max_az_idx);
estimated_elevation = fine_el_grid(max_el_idx);

% 尝试搜索第二个峰值（可能对应反射路径）
mask = fine_spectrum;
% 将最大峰值附近的点置零，寻找次高峰
mask_radius = 3;  % 掩码半径
for i = max(1, max_az_idx-mask_radius):min(fine_az_steps, max_az_idx+mask_radius)
    for j = max(1, max_el_idx-mask_radius):min(fine_el_steps, max_el_idx+mask_radius)
        mask(i, j) = 0;
    end
end
[second_peak_value, second_max_idx] = max(mask(:));
[second_max_az_idx, second_max_el_idx] = ind2sub(size(mask), second_max_idx);
second_az = fine_az_grid(second_max_az_idx);
second_el = fine_el_grid(second_max_el_idx);

% 打印主要和次要峰值
fprintf('MUSIC主峰: 方位角=%.2f°, 俯仰角=%.2f°, 峰值=%.2e\n', ...
        estimated_azimuth, estimated_elevation, peak_value);
if second_peak_value > 0.5 * peak_value
    fprintf('MUSIC次峰: 方位角=%.2f°, 俯仰角=%.2f°, 峰值=%.2e (%.1f%%的主峰)\n', ...
            second_az, second_el, second_peak_value, second_peak_value/peak_value*100);
end

% 调整俯仰角的符号，使其与仿真坐标系一致
% 在该系统中，正z对应正俯仰角
estimated_elevation = -estimated_elevation;

% 更新下一次调用的先验信息
params.est.prior_azimuth = estimated_azimuth;
params.est.prior_elevation = estimated_elevation;
params.est.prior_az_std = max(5, az_std * 0.8);  % 逐渐减小标准差，对结果更有信心
params.est.prior_el_std = max(5, el_std * 0.8);

end

function a = steering_vector(az, el, rx_elements_pos, params)
    % 计算给定方位角和俯仰角的导向矢量
    lambda = params.c / params.fc;  % 波长
    
    % 计算单位方向向量
    az_rad = deg2rad(az);
    el_rad = deg2rad(el);
    k_vec = 2*pi/lambda * [cosd(el)*cosd(az); cosd(el)*sind(az); sind(el)];
    
    % 计算每个阵元的空间相位
    % rx_elements_pos是Nx3矩阵，每行表示一个阵元的[x,y,z]坐标
    % 我们需要计算每个阵元的相位: exp(j * k·r)
    phase = zeros(size(rx_elements_pos, 1), 1);
    for i = 1:size(rx_elements_pos, 1)
        % 计算k·r (点积)
        phase(i) = k_vec(1) * rx_elements_pos(i, 1) + ...
                   k_vec(2) * rx_elements_pos(i, 2) + ...
                   k_vec(3) * rx_elements_pos(i, 3);
    end
    a = exp(1j * phase);
    
    % 归一化导向矢量
    a = a / norm(a);
end 