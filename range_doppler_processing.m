function [range_doppler, range_axis, doppler_axis] = range_doppler_processing(rx_signal, params)
%RANGE_DOPPLER_PROCESSING 对接收信号进行距离-多普勒处理
%   rx_signal: 接收信号矩阵 [采样点数 x chirp数 x 接收阵元数]
%   params: 系统参数结构体
%   range_doppler: 距离-多普勒谱 [距离点数 x 多普勒点数]
%   range_axis: 距离轴 (m)
%   doppler_axis: 多普勒轴 (Hz)

% 提取参数
c = params.c;                   % 光速
lambda = params.lambda;         % 波长
fs = params.fmcw.fs;            % 采样率
T = params.fmcw.T;              % 扫频时间
B = params.fmcw.B;              % 带宽
num_chirps = params.fmcw.num_chirps;  % chirp数量
nfft_range = params.rd.nfft_range;    % 距离维FFT点数
nfft_doppler = params.rd.nfft_doppler;  % 多普勒维FFT点数
window_range = params.rd.window_range;    % 距离维窗函数
window_doppler = params.rd.window_doppler;  % 多普勒维窗函数
sweep_rate = params.fmcw.mu;   % 调频率 (使用sweep_rate代替mu)

% 获取接收信号维度
[num_samples, ~, num_rx] = size(rx_signal);

% 确保nfft_range适当大小，不要过度扩展
if nfft_range > 4 * num_samples
    warning('nfft_range过大，可能导致距离估计失真。调整为2倍采样点数');
    nfft_range = 2 * num_samples;
    if mod(nfft_range, 2) ~= 0
        nfft_range = nfft_range + 1;  % 确保是偶数
    end
end

% 每个距离bin在FFT后的点数
num_range_bins = nfft_range/2;

% 对所有接收天线进行处理并平均
range_doppler_all = zeros(num_range_bins, nfft_doppler, num_rx);

for rx_idx = 1:num_rx
    % 获取当前接收天线的信号
    current_rx_signal = rx_signal(:, :, rx_idx);
    
    % 1. 距离FFT (快时间FFT)
        % 应用窗函数
    if strcmp(window_range, 'hamming')
        window_r = hamming(num_samples);
    elseif strcmp(window_range, 'hanning')
        window_r = hanning(num_samples);
    else
        window_r = ones(num_samples, 1);  % 矩形窗/无窗
    end
    
    % 应用窗函数到每个chirp
    windowed_signal = current_rx_signal .* repmat(window_r, 1, num_chirps);
    
    % 对每个chirp做FFT (沿快时间维度)
    range_fft = fft(windowed_signal, nfft_range, 1);
    
    % 只保留前一半频率点 (负频率是镜像)
    range_fft = range_fft(1:num_range_bins, :);
    
    % 2. 多普勒FFT (慢时间FFT)
    % 应用窗函数
    if strcmp(window_doppler, 'hamming')
        window_d = hamming(num_chirps);
    elseif strcmp(window_doppler, 'hanning')
        window_d = hanning(num_chirps);
    else
        window_d = ones(num_chirps, 1);  % 矩形窗/无窗
    end
    
    % 应用窗函数到每个距离bin
    windowed_range = range_fft .* repmat(window_d', num_range_bins, 1);

    % 对每个距离bin做FFT (沿慢时间维度)
    range_doppler_tmp = fftshift(fft(windowed_range, nfft_doppler, 2), 2);
    
    % 保存当前接收天线的处理结果 - 确保维度匹配
    if size(range_doppler_tmp, 1) == num_range_bins && size(range_doppler_tmp, 2) == nfft_doppler
        range_doppler_all(:, :, rx_idx) = abs(range_doppler_tmp);
    else
        % 如果维度不匹配，输出错误信息并调整大小
        fprintf('维度不匹配: range_doppler_tmp尺寸 %dx%d, 期望 %dx%d\n', ...
            size(range_doppler_tmp, 1), size(range_doppler_tmp, 2), num_range_bins, nfft_doppler);
        
        % 创建正确尺寸的临时矩阵
        temp = zeros(num_range_bins, nfft_doppler);
        
        % 复制可用数据
        copy_rows = min(size(range_doppler_tmp, 1), num_range_bins);
        copy_cols = min(size(range_doppler_tmp, 2), nfft_doppler);
        temp(1:copy_rows, 1:copy_cols) = abs(range_doppler_tmp(1:copy_rows, 1:copy_cols));
        
        range_doppler_all(:, :, rx_idx) = temp;
    end
end

% 对所有接收天线的结果求平均 (非相干积累)
range_doppler = mean(range_doppler_all, 3);

% 计算距离轴
% FMCW雷达中距离与拍频的正确关系: R = (f_beat * c) / (2 * μ)
% 频率分辨率 = fs/nfft_range
% 距离分辨率 = c / (2 * B)
range_res = c / (2 * B);

% 频率轴
freq_axis = (0:num_range_bins-1) * (fs/nfft_range);

% 正确的距离轴计算，基于频率到距离的转换
range_axis = freq_axis * c / (2 * sweep_rate);

% 计算理论最大探测距离（基于采样率）
max_unambiguous_range = (fs/2) * c / (2 * sweep_rate);
fprintf('距离处理参数: 采样点数=%d, FFT点数=%d\n', num_samples, nfft_range);
fprintf('距离分辨率=%.2f m, 最大无模糊距离=%.2f m\n', range_res, max_unambiguous_range);

% 计算多普勒轴
% 多普勒分辨率 = 2 / (lambda * N * T)，其中N是chirp数
doppler_res = 2 / (lambda * num_chirps * T);
doppler_axis = (-nfft_doppler/2:nfft_doppler/2-1) * doppler_res;

end 