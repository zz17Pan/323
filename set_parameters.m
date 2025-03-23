function params = set_parameters()
%SET_PARAMETERS 配置系统全部参数
%   返回包含所有系统参数的结构体

%% 基本参数
params.c = 3e8;                  % 光速(m/s)
params.fc = 100e9;               % 载波频率 (100 GHz)
params.lambda = params.c / params.fc;  % 波长 (m)

%% FMCW信号参数
params.fmcw.T = 1e-3;            % 扫频时间 (s)
params.fmcw.B = 4e9;             % 带宽 (增加到4 GHz，提高距离分辨率)
params.fmcw.mu = params.fmcw.B / params.fmcw.T;  % 调频率
params.fmcw.fs = 21e6;           % 采样率(20 MHz，设置为合理值以降低计算负担)
% 确保采样时间覆盖整个扫频时间
params.fmcw.Ns = ceil(params.fmcw.T * params.fmcw.fs);  % 计算至少需要的采样点数
fprintf('更新采样点数以覆盖扫频时间 Ns = %d\n', params.fmcw.Ns);

% 计算理论最大距离和分辨率
max_theoretical_range = params.c * params.fmcw.fs / (2 * params.fmcw.mu * 2); % 理论最大探测距离(单位：米)
range_resolution = params.c / (2 * params.fmcw.B); % 距离分辨率(单位：米)
fprintf('FMCW理论参数 - 最大探测距离 %.2f m, 距离分辨率 %.2f m\n', max_theoretical_range, range_resolution);

% 确保采样点数是偶数，便于FFT处理
if mod(params.fmcw.Ns, 2) ~= 0
    params.fmcw.Ns = params.fmcw.Ns + 1;
end

params.fmcw.num_chirps = 128;    % 每帧chirp数
params.fmcw.A = 1;               % 发射信号幅度

%% 发射阵列参数
params.tx.array_size = [4, 4];   % 发射阵列大小 [行 列]
params.tx.spacing = params.lambda / 2;  % 阵元间距 (半波长)
params.tx.pos = [0, 0, 0];       % 发射阵列位置 (坐标原点)

%% 接收阵列参数
params.rx.array_size = [4, 4];   % 接收阵列大小 [行 列]
params.rx.spacing = params.lambda / 2;  % 阵元间距 (半波长)
% 修改初始位置，避开特殊角度0度和90度
params.rx.init_pos = [10, 15, 100];  % 接收阵列初始位置 (避免y=0的特殊情况)
params.rx.velocity = [1, 0.5, 0.2];  % 接收阵列速度 (m/s) [vx, vy, vz]

%% 信道模型参数
params.channel.add_noise = true;  % 是否添加噪声
params.channel.snr = 20;         % 信噪比(dB)
params.channel.num_reflectors = 1;  % 反射体数(多径)，减少干扰
params.channel.reflection_coef = 0.1;  % 反射系数，降低反射强度

%% 距离-多普勒处理参数
params.rd.window_range = 'hamming';   % 距离维窗函数
params.rd.window_doppler = 'hamming'; % 多普勒维窗函数

% 确保FFT点数是采样点数的合理倍数，避免过大导致距离失真
% 使用2倍采样点数是更合理的选择
params.rd.nfft_range = 2 * params.fmcw.Ns;  % 距离维FFT点数修改为2倍采样点数
if mod(params.rd.nfft_range, 2) ~= 0
    params.rd.nfft_range = params.rd.nfft_range + 1;  % 确保是偶数
end
params.rd.nfft_doppler = 2^8;        % 多普勒维FFT点数

%% CFAR参数
params.cfar.guard_cells = [4, 4];     % 保护单元 [距离, 多普勒]
params.cfar.training_cells = [8, 8];  % 训练单元 [距离, 多普勒]
params.cfar.pfa = 1e-5;              % 虚警率，提高检测灵敏度
params.cfar.threshold_factor = 2.5;   % CFAR检测阈值系数，值越大检测越保守

%% MUSIC参数
params.music.num_sources = 1;       % 信源数
params.music.az_range = [-90, 90];  % 方位角搜索范围(-90度, 90度)，扩大搜索范围
params.music.el_range = [-90, 90];  % 俯仰角搜索范围(-90度, 90度)，扩大搜索范围
params.music.az_resolution = 2;     % 方位角分辨率 (度)，降低分辨率以减少计算量
params.music.el_resolution = 2;     % 俯仰角分辨率 (度)，降低分辨率以减少计算量

%% OMP参数
params.omp.max_iter = 20;           % 最大迭代次数
params.omp.residual_tol = 0.01;     % 残差阈值
params.omp.range_grid_size = 0.5;   % 距离网格大小 (m)
params.omp.angle_grid_size = 2;     % 角度网格大小 (度)，增大网格尺寸以减少计算量
params.omp.prior_scale = 5;         % 先验分布缩放因子 (搜索范围±5σ)，扩大搜索范围

%% 估计误差参数
params.est.range_var = 1^2;         % 距离估计方差 (m^2)，增大容错范围
params.est.azimuth_var = 2^2;       % 方位角估计方差(度^2)，增大容错范围
params.est.elevation_var = 2^2;     % 俯仰角估计方差(度^2)，增大容错范围

%% 卡尔曼滤波参数
params.kf.dt = 0.1;                 % 采样间隔 (s)
params.kf.q_pos = 0.5^2;            % 位置过程噪声方差，增大以适应更快的收敛
params.kf.q_vel = 0.05^2;           % 速度过程噪声方差，增大以适应更快的收敛
params.kf.r_range = 2^2;            % 距离测量噪声方差，增大以适应初始误差
params.kf.r_azimuth = 5^2;          % 方位角测量噪声方差，增大以适应初始误差
params.kf.r_elevation = 5^2;        % 俯仰角测量噪声方差，增大以适应初始误差

%% 仿真参数
params.sim.num_frames = 20;         % 帧数，减少帧数以加快测试
params.sim.frame_interval = 0.1;    % 帧间间隔(s)

%% 可视化参数
params.viz.update_interval = 2;     % 可视化更新间隔(秒)

end 
