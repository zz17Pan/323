function rx_signal = simulate_propagation(tx_signal, tx_array, rx_array, params)
%SIMULATE_PROPAGATION 模拟信号从发射端传播到接收端的过程
%   tx_signal: 发射信号矩阵 [采样点数 x chirp数]
%   tx_array: 发射阵列结构体
%   rx_array: 接收阵列结构体
%   params: 系统参数结构体
%   rx_signal: 接收信号矩阵 [采样点数 x chirp数 x 接收阵元数]

% 提取参数
c = params.c;               % 光速
lambda = params.lambda;     % 波长
fc = params.fc;             % 载波频率
T = params.fmcw.T;          % 扫频时间
fs = params.fmcw.fs;        % 采样率
Ns = params.fmcw.Ns;        % 每个chirp的采样点数
num_chirps = params.fmcw.num_chirps;  % chirp数量
sweep_rate = params.fmcw.mu;  % 调频率

% 计算接收阵列相对于发射阵列的球坐标参数
rx_center = rx_array.pos;  % 接收阵列中心

% 计算水平距离和总距离
horizontal_distance = sqrt(rx_center(1)^2 + rx_center(2)^2);
r = sqrt(rx_center(1)^2 + rx_center(2)^2 + rx_center(3)^2);

% 计算方位角 (水平面内从x轴顺时针方向的角度)
azimuth = atan2d(rx_center(2), rx_center(1));

% 计算俯仰角 (从水平面到目标的角度，向上为正)
elevation = atan2d(rx_center(3), horizontal_distance);

% 打印详细的角度信息以便调试
fprintf('模拟信号传播 - 接收中心: [%.2f, %.2f, %.2f], 距离: %.2f m\n', rx_center(1), rx_center(2), rx_center(3), r);
fprintf('计算的方位角: %.2f°, 俯仰角: %.2f°\n', azimuth, elevation);

% 计算径向速度 (使用连续两帧的位置估计)
if isfield(rx_array, 'prev_pos')
    delta_t = params.sim.frame_interval;
    velocity = (rx_center - rx_array.prev_pos) / delta_t;
    % 径向速度是速度在发射端到接收端连线方向上的投影
    direction = rx_center / norm(rx_center);  % 单位方向向量
    radial_velocity = dot(velocity, direction);
else
    % 初始帧没有前一个位置，假设径向速度为0
    radial_velocity = 0;
end
rx_array.prev_pos = rx_center;  % 保存当前位置用于下一帧计算

% 计算多普勒频移
doppler_shift = 2 * radial_velocity * fc / c;

% 生成HSPM信道矩阵 (球面波-平面波混合模型)
% 这里简化为直接计算每个接收阵元对应的时延和多普勒
num_rx_elements = rx_array.num_elements;
rx_signal = zeros(Ns, num_chirps, num_rx_elements);

% 时间向量
t = (0:Ns-1)' / fs;  % 列向量

% 预计算一些值以加速处理
chirp_times = (0:num_chirps-1) * T;  % 每个chirp的时间偏移

% 将角度转为弧度进行计算
az_rad = deg2rad(azimuth);
el_rad = deg2rad(elevation);

% 计算波数向量 (与MUSIC中保持一致)
% x = r·cos(el)·cos(az)
% y = r·cos(el)·sin(az)
% z = r·sin(el)
k_vector = 2*pi/lambda * [cos(el_rad)*cos(az_rad); 
                          cos(el_rad)*sin(az_rad); 
                          sin(el_rad)];

% LoS路径 (主路径)
for rx_idx = 1:num_rx_elements
    % 计算发射端到当前接收阵元的距离
    rx_element_pos = rx_array.elements_pos(rx_idx, :);
    distance = norm(rx_element_pos);
    
    % 计算该接收元素的方向向量
    direction = rx_element_pos / norm(rx_element_pos);  % 单位方向向量
    
    % 计算时延
    tau = distance / c;
    
    % 时延是否在合理范围内 (避免处理不必要的超远距离)
    if tau > T
        continue;  % 时延超过一个chirp的持续时间，跳过处理
    end
    
    % 计算相位偏移 (和接收阵元位置有关)
    pos_vector = rx_element_pos'; % 转为列向量
    phase_shift = k_vector' * pos_vector; % 空间相位
    
    % 计算拍频
    beat_freq = sweep_rate * tau;  % 拍频频率
    
    % 计算时移后的基带信号
    time_shift = t - tau;
    
    % 只处理时延在采样窗口内的部分
    valid_samples = time_shift >= 0 & time_shift < T;
    
    % 如果没有有效样本，跳过
    if ~any(valid_samples)
        continue;
    end
    
    % 计算延时信号 (基带信号)
    delayed_signal = zeros(Ns, 1);
    delayed_signal(valid_samples) = exp(1j * 2 * pi * (0.5 * sweep_rate * time_shift(valid_samples).^2));
    
    % 计算多普勒相位因子
    doppler_phase_factor = exp(1j * 2 * pi * doppler_shift * t);
    
    % 空间相位因子
    spatial_phase_factor = exp(1j * phase_shift);
    
    % 对每个chirp计算接收信号
    for chirp_idx = 1:num_chirps
        % 当前chirp的时间偏移
        chirp_time = chirp_times(chirp_idx);
        
        % 计算当前chirp的总多普勒相位
        chirp_doppler_phase = exp(1j * 2 * pi * doppler_shift * chirp_time);
        
        % 组合各种相位因子
        rx_signal(:, chirp_idx, rx_idx) = delayed_signal .* doppler_phase_factor * (chirp_doppler_phase * spatial_phase_factor);
    end
end

% 加入NLoS路径 (反射) - 优化性能
if params.channel.num_reflectors > 0
    % 这里简化为随机生成几个反射路径
    for reflector_idx = 1:params.channel.num_reflectors
        % 随机生成反射体位置，但限制在更合理的范围内
        max_reflector_range = min(50, norm(rx_array.pos)/4);  % 反射体最大距离
        reflector_pos = [(rand-0.5)*max_reflector_range, (rand-0.5)*max_reflector_range, max_reflector_range/2 + (rand-0.5)*max_reflector_range/4];
        
        % 反射路径的衰减系数
        attenuation = params.channel.reflection_coef;
        
        % 为该反射体处理所有接收阵元
        for rx_idx = 1:num_rx_elements
            % 计算反射路径的总距离: 发射端->反射体->接收阵元
            rx_element_pos = rx_array.elements_pos(rx_idx, :);
            distance_tx_to_reflector = norm(reflector_pos);
            distance_reflector_to_rx = norm(reflector_pos - rx_element_pos);
            total_distance = distance_tx_to_reflector + distance_reflector_to_rx;
            
            % 计算反射路径时延
            tau_reflect = total_distance / c;
            
            % 跳过超出合理范围的时延
            if tau_reflect > T
                continue;
            end
            
            % 计算反射路径的相位偏移 (简化)
            phase_shift_reflect = 2*pi * rand;  % 随机相位偏移
            
            % 计算时移后的基带信号
            time_shift = t - tau_reflect;
            
            % 只处理时延在采样窗口内的部分
            valid_samples = time_shift >= 0 & time_shift < T;
            
            % 如果没有有效样本，跳过
            if ~any(valid_samples)
                continue;
            end
            
            % 计算延时信号 (基带信号)
            delayed_signal = zeros(Ns, 1);
            delayed_signal(valid_samples) = exp(1j * 2 * pi * (0.5 * sweep_rate * time_shift(valid_samples).^2));
            
            % 考虑反射体可能的移动带来的多普勒效应 (简化为与主路径相同)
            % 对每个chirp计算接收信号
            for chirp_idx = 1:num_chirps
                % 当前chirp的时间偏移
                chirp_time = chirp_times(chirp_idx);
                
                % 计算当前chirp的多普勒相位
                chirp_doppler_phase = exp(1j * 2 * pi * doppler_shift * chirp_time);
                
                % 计算反射信号
                reflection = attenuation * delayed_signal * chirp_doppler_phase * exp(1j * phase_shift_reflect);
                
                % 将反射信号加到接收信号上
                rx_signal(:, chirp_idx, rx_idx) = rx_signal(:, chirp_idx, rx_idx) + reflection;
            end
        end
    end
end

% 添加噪声
if params.channel.add_noise
    % 计算信号功率
    signal_power = mean(abs(rx_signal(:)).^2);
    
    % 如果信号功率接近零，设置一个最小值以避免除零错误
    if signal_power < 1e-10
        signal_power = 1e-10;
    end
    
    % 根据SNR计算噪声功率
    snr_linear = 10^(params.channel.snr/10);
    noise_power = signal_power / snr_linear;
    
    % 生成复高斯白噪声
    noise = sqrt(noise_power/2) * (randn(size(rx_signal)) + 1j*randn(size(rx_signal)));
    
    % 将噪声加入接收信号
    rx_signal = rx_signal + noise;
end

end 