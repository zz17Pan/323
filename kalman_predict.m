function kf = kalman_predict(kf, params)
%KALMAN_PREDICT 执行卡尔曼滤波器的预测步骤
%   kf: 卡尔曼滤波器状态结构体
%   params: 系统参数结构体
%   更新后的卡尔曼滤波器状态

% 状态预测
kf.x = kf.A * kf.x;

% 协方差预测
kf.P = kf.A * kf.P * kf.A' + kf.Q;

end 