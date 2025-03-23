function tx_signal = generate_fmcw(params)
%GENERATE_FMCW 生成FMCW发射信号
%   params: 系统参数结构体
%   tx_signal: 发射信号矩阵 [采样点数 x chirp数]

% 提取参数
T = params.fmcw.T;         % 扫频时间
fs = params.fmcw.fs;       % 采样率
Ns = params.fmcw.Ns;       % 每个chirp的采样点数
num_chirps = params.fmcw.num_chirps;  % chirp数量
fc = params.fc;            % 载波频率
sweep_rate = params.fmcw.mu;       % 调频率 (B/T)
A = params.fmcw.A;         % 信号幅度

% 时间向量 (每个chirp的)
t = (0:Ns-1)' / fs;  % 列向量
fprintf('FMCW信号生成 - 采样点数: %d, 时间范围: [0, %.3f] ms\n', Ns, max(t)*1000);

% 生成基带FMCW信号
% s_tx(t) = A * exp(j*2*pi*(fc*t + 0.5*mu*t^2))
% 不考虑载波fc (因为基带信号处理)
baseband_signal = A * exp(1j * 2 * pi * (0.5 * sweep_rate * t.^2));

% 检查信号是否覆盖了足够的时间
% 使用1%的容差，避免浮点数比较精度问题
if max(t)*1.01 < T
    warning('采样时间 (%.3f ms) 小于扫频时间 (%.3f ms)，可能导致距离估计不准确。', max(t)*1000, T*1000);
else
    fprintf('采样时间充分覆盖了扫频时间\n');
end

% 预分配内存
tx_signal = zeros(Ns, num_chirps);

% 对每个chirp使用相同的基带信号
for i = 1:num_chirps
    tx_signal(:, i) = baseband_signal;
end

end 