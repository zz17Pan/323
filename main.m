%% 太赫兹波束对准感知部分主函数
% 作者：
% 日期：

close all;
clear;
clc;

%% 参数设置
params = set_parameters();

% 打印关键参数
fprintf('FMCW参数: T=%.3f ms, B=%.2f GHz, fs=%.2f MHz, 采样点数=%d\n', ...
    params.fmcw.T*1e3, params.fmcw.B/1e9, params.fmcw.fs/1e6, params.fmcw.Ns);
fprintf('阵列配置: 发射端 %dx%d, 接收端 %dx%d\n', ...
    params.tx.array_size(1), params.tx.array_size(2), ...
    params.rx.array_size(1), params.rx.array_size(2));

%% 初始化发射接收阵列
try
    [tx_array, rx_array] = init_arrays(params);
    fprintf('阵列初始化成功: 发射端 %d 个阵元, 接收端 %d 个阵元\n', ...
        tx_array.num_elements, rx_array.num_elements);
catch ME
    fprintf('阵列初始化错误: %s\n', ME.message);
    rethrow(ME);
end

%% 主循环 - 模拟接收端移动并进行感知
num_frames = params.sim.num_frames;
results = struct('time', cell(1, num_frames), 'true', cell(1, num_frames), 'est', cell(1, num_frames));

% 卡尔曼滤波器初始化
kf = init_kalman_filter(params);

% 接收端初始位置和速度
rx_pos = params.rx.init_pos;
rx_vel = params.rx.velocity;

% 设置计时器以估计剩余时间
total_timer = tic;
frame_times = zeros(1, num_frames);

for frame_idx = 1:num_frames
    frame_timer = tic;
    fprintf('\n处理帧 %d/%d...\n', frame_idx, num_frames);
    
    % 更新接收端位置
    time = (frame_idx-1) * params.sim.frame_interval;
    rx_pos = params.rx.init_pos + rx_vel * time;
    
    % 计算接收端相对于发射端的球坐标
    % 计算水平距离和总距离
    horizontal_distance = sqrt(rx_pos(1)^2 + rx_pos(2)^2);
    range = sqrt(rx_pos(1)^2 + rx_pos(2)^2 + rx_pos(3)^2);
    
    % 计算方位角 (水平面内从x轴顺时针方向的角度)
    azimuth = atan2d(rx_pos(2), rx_pos(1));
    
    % 计算俯仰角 (从水平面到目标的角度，向上为正)
    elevation = atan2d(rx_pos(3), horizontal_distance);
    
    fprintf('真实位置: 距离=%.2f m, 方位角=%.2f°, 俯仰角=%.2f°\n', range, azimuth, elevation);
    fprintf('接收端坐标: x=%.2f m, y=%.2f m, z=%.2f m\n', rx_pos(1), rx_pos(2), rx_pos(3));
    
    % 更新接收端阵列位置
    rx_array = update_rx_array(rx_array, rx_pos);
    
    % 添加异常处理，防止单个帧的错误导致整个程序崩溃
    try
        % 预测阶段 - 使用上一帧的结果预测当前帧的位置
        % 注意：第一帧没有先验信息，使用初始设置
        if frame_idx == 1
            % 卡尔曼滤波预测，没有先验时的特殊处理
            kf = kalman_predict(kf, params);
            
            % 第一帧的先验信息从卡尔曼滤波获取
            prior_info.range = kf.x(1);
            prior_info.azimuth = kf.x(3);
            prior_info.elevation = kf.x(5);
            
            % 并扩大初始搜索协方差，以适应初始状态的不确定性
            prior_cov = diag([params.est.range_var*4, params.est.azimuth_var*4, params.est.elevation_var*4]);
        else
            % 正常的卡尔曼预测
            kf = kalman_predict(kf, params);
            
            % 使用卡尔曼滤波预测作为先验信息
            prior_info.range = kf.x(1);
            prior_info.azimuth = kf.x(3);
            prior_info.elevation = kf.x(5);
            
            % 根据卡尔曼滤波协方差设置合适的搜索范围
            prior_cov = diag([kf.P(1,1), kf.P(3,3), kf.P(5,5)]);
        end
        fprintf('卡尔曼滤波预测: 距离=%.2f m, 方位角=%.2f°, 俯仰角=%.2f°\n', prior_info.range, prior_info.azimuth, prior_info.elevation);
        
        % 为OMP算法准备额外参数
        params_omp = params;
        params_omp.omp.max_iter = 3;  % 减少迭代次数以提高速度
        
        % 生成发射信号
        tx_signal = generate_fmcw(params);
        
        % 模拟传播并得到接收信号
        rx_signal = simulate_propagation(tx_signal, tx_array, rx_array, params);
        
        % 距离多普勒处理
        [rd_cube, range_axis, velocity_axis] = range_doppler_processing(rx_signal, params);
        
        % CFAR检测，获取距离和速度
        [detected_range, detected_velocity] = cfar_detection(rd_cube, range_axis, velocity_axis, params, prior_info.range);
        
        % 使用CFAR结果作为先验距离
        if ~isempty(detected_range)
            % 将检测到的单个距离值赋给先验
            prior_info.range = detected_range;
            
            % 选取与当前处理峰值对应的采样点范围
            peak_idx = round(prior_info.range * 2 * params.fmcw.mu / params.c * params.fmcw.fs);
            peak_idx = max(1, min(peak_idx, size(rx_signal, 1)-1));  % 确保索引有效
            
            fprintf('CFAR检测: 距离=%.2f m, 速度=%.2f m/s\n', detected_range, detected_velocity);
        else
            fprintf('CFAR检测: 未检测到有效目标\n');
        end
        
        % OMP稀疏重建 - 结合优化的先验信息
        [est_range, est_azimuth, est_elevation] = prior_guided_omp(rx_signal, tx_array, rx_array, prior_info, prior_cov, params_omp);
        fprintf('OMP估计结果: 距离=%.2f m, 方位角=%.2f°, 俯仰角=%.2f°\n', est_range, est_azimuth, est_elevation);
        
        % 卡尔曼滤波更新
        z = [est_range; est_azimuth; est_elevation];
        kf = kalman_update(kf, z);
        fprintf('卡尔曼滤波结果: 距离=%.2f m, 方位角=%.2f°, 俯仰角=%.2f°\n', kf.x(1), kf.x(3), kf.x(5));
        
        % 存储结果
        results(frame_idx).time = time;
        results(frame_idx).true = [range, azimuth, elevation];
        results(frame_idx).est = [kf.x(1), kf.x(3), kf.x(5)]; % 滤波后的估计值
        
        % 计算并显示当前帧误差
        range_error = kf.x(1) - range;
        az_error = kf.x(3) - azimuth;
        el_error = kf.x(5) - elevation;
        
        fprintf('\n当前帧估计误差:\n');
        fprintf('距离: 真实=%.2f m, 估计=%.2f m, 误差=%.2f m, 比例=%.2f\n', ...
                range, kf.x(1), range_error, range_error/range);
        fprintf('方位角: 真实=%.2f°, 估计=%.2f°, 误差=%.2f°\n', ...
                azimuth, kf.x(3), az_error);
        fprintf('俯仰角: 真实=%.2f°, 估计=%.2f°, 误差=%.2f°\n', ...
                elevation, kf.x(5), el_error);
            
        % 可视化当前结果
        if mod(frame_idx, params.viz.update_interval) == 0
            visualize_tracking(results(1:frame_idx), params);
        end
        
    catch ME
        % 捕获错误并打印详细信息
        fprintf('错误发生在帧 %d: %s\n', frame_idx, ME.message);
        fprintf('堆栈追踪:\n');
        for k = 1:length(ME.stack)
            fprintf('  文件: %s, 行: %d, 函数: %s\n', ...
                ME.stack(k).file, ME.stack(k).line, ME.stack(k).name);
        end
        fprintf('已跳过问题帧，继续处理\n');
        
        % 保持卡尔曼滤波器状态不变，避免对后续帧产生影响
        % 或者可以使用预测值作为当前帧的估计结果
        results(frame_idx).time = time;
        results(frame_idx).true = [range, azimuth, elevation];
        
        % 如果是第一帧出错，则使用真实值进行初始化
        if frame_idx == 1
            kf.x = [range; 0; azimuth; 0; elevation; 0];
            results(frame_idx).est = [range, azimuth, elevation];
        else
            % 使用预测值作为估计结果
            results(frame_idx).est = [kf.x(1), kf.x(3), kf.x(5)];
        end
    end
    
    % 记录帧处理时间
    frame_times(frame_idx) = toc(frame_timer);
    mean_frame_time = mean(frame_times(1:frame_idx));
    est_remaining_time = mean_frame_time * (num_frames - frame_idx);
    
    fprintf('帧处理时间: %.2f 秒, 估计剩余时间: %.1f 分钟\n', ...
        frame_times(frame_idx), est_remaining_time/60);
end

total_time = toc(total_timer);
fprintf('\n处理完成。总耗时: %.1f 分钟\n', total_time/60);

%% 结果评估
try
    evaluate_performance(results, params);
catch ME
    fprintf('结果评估错误: %s\n', ME.message);
end