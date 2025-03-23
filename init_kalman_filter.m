function kf = init_kalman_filter(params)
%INIT_KALMAN_FILTER 初始化卡尔曼滤波器
%   params: 系统参数结构体
%   kf: 卡尔曼滤波器结构体

% 状态向量 x = [r; vr; az; vaz; el; vel]
% r: 距离
% vr: 径向速度
% az: 方位角
% vaz: 方位角速度
% el: 俯仰角
% vel: 俯仰角速度

% 从参数计算初始状态
initial_pos = params.rx.init_pos;
initial_vel = params.rx.velocity;

% 计算初始距离和角度
r0 = norm(initial_pos);
horizontal_dist = sqrt(initial_pos(1)^2 + initial_pos(2)^2);
az0 = atan2d(initial_pos(2), initial_pos(1));
el0 = atan2d(initial_pos(3), horizontal_dist);

% 计算径向速度
r_dir = initial_pos / r0;  % 径向单位向量
vr0 = dot(initial_vel, r_dir);  % 径向速度分量

% 计算角速度（使用更准确的方法）
if norm(initial_vel) > 0
    % 使用叉积计算切向速度，然后计算角速度
    cross_product = cross(initial_pos, initial_vel);
    tangential_vel = norm(cross_product) / r0;
    
    % 方位角速度（水平面内）
    vaz0 = tangential_vel / (horizontal_dist + 1e-10) * sign(cross_product(3));
    
    % 俯仰角速度（垂直面内）
    horizontal_vel = sqrt(initial_vel(1)^2 + initial_vel(2)^2);
    vel0 = (initial_vel(3) * horizontal_dist - initial_pos(3) * horizontal_vel) / (r0^2 + 1e-10);
else
    vaz0 = 0;
    vel0 = 0;
end

% 初始状态向量
x0 = [r0; vr0; az0; vaz0; el0; vel0];

% 状态转移矩阵 - 使用更精确的非线性模型
dt = params.kf.dt;  % 时间步长

% 调整状态转移矩阵
A = [1, dt, 0, 0, 0, 0;    % r' = r + vr*dt
     0, 0.95, 0, 0, 0, 0;     % vr' = 0.95*vr（轻微阻尼）
     0, 0, 1, dt, 0, 0;    % az' = az + vaz*dt
     0, 0, 0, 0.95, 0, 0;  % vaz' = 0.95*vaz（轻微阻尼）
     0, 0, 0, 0, 1, dt;    % el' = el + vel*dt
     0, 0, 0, 0, 0, 0.95]; % vel' = 0.95*vel（轻微阻尼）

% 观测矩阵（目前只直接观测距离、方位角和俯仰角）
H = [1, 0, 0, 0, 0, 0;     % 观测距离
     0, 0, 1, 0, 0, 0;     % 观测方位角
     0, 0, 0, 0, 1, 0];    % 观测俯仰角

% 重新调整过程噪声协方差
q_pos = params.kf.q_pos * 2.0;    % 位置过程噪声方差增大2倍
q_vel = params.kf.q_vel * 3.0;    % 速度过程噪声方差增大3倍

% 调整不同状态维度的过程噪声
Q = diag([
    q_pos*5,    % 距离噪声适度增大
    q_vel*3,    % 径向速度噪声增大
    q_pos*1.5,    % 方位角噪声适度增大
    q_vel*3,    % 方位角速度噪声增大
    q_pos*1.5,    % 俯仰角噪声适度增大
    q_vel*3     % 俯仰角速度噪声增大
]);

% 测量噪声协方差 - 调整为更合理的值
r_range = params.kf.r_range * 0.8;       % 距离测量噪声降低20%
r_azimuth = params.kf.r_azimuth * 0.8;   % 方位角测量噪声降低20%
r_elevation = params.kf.r_elevation * 0.8; % 俯仰角测量噪声降低20%

% 基于实际测量表现调整测量噪声
R = diag([
    r_range,     % 距离噪声
    r_azimuth,     % 方位角噪声
    r_elevation    % 俯仰角噪声
]);

% 初始状态估计协方差 - 反映初始估计的不确定性
% 考虑不同状态分量的尺度差异
P0 = eye(6);
P0(1,1) = r_range * 2;    % 距离初始不确定性
P0(2,2) = q_vel * 2;      % 径向速度初始不确定性
P0(3,3) = r_azimuth * 2;  % 方位角初始不确定性
P0(4,4) = q_vel;          % 方位角速度初始不确定性
P0(5,5) = r_elevation * 2;% 俯仰角初始不确定性
P0(6,6) = q_vel;          % 俯仰角速度初始不确定性

% 创建卡尔曼滤波器结构体
kf = struct('x', x0, 'P', P0, 'A', A, 'H', H, 'Q', Q, 'R', R);

% 打印初始状态
fprintf('卡尔曼滤波器初始化: 初始状态 [r=%.2f, vr=%.2f, az=%.2f, vaz=%.2f, el=%.2f, vel=%.2f]\n', ...
    x0(1), x0(2), x0(3), x0(4), x0(5), x0(6));

end 