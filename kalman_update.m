function kf = kalman_update(kf, z)
%KALMAN_UPDATE 卡尔曼滤波器更新步骤
%   kf: 卡尔曼滤波器结构体
%   z: 观测向量 [距离; 方位角; 俯仰角]
%   返回更新后的卡尔曼滤波器结构体

% 计算观测预测
z_pred = kf.H * kf.x;

% 计算观测残差（通常称为创新）
y = z - z_pred;

% 角度归一化（解决角度绕行问题）
y(2) = wrapToPi(deg2rad(y(2))) * 180/pi;  % 方位角归一化到[-180, 180]
y(3) = wrapToPi(deg2rad(y(3))) * 180/pi;  % 俯仰角归一化到[-180, 180]

% 计算创新协方差
S = kf.H * kf.P * kf.H' + kf.R;

% 计算卡尔曼增益
K = kf.P * kf.H' / S;

% 检测异常值（过滤不合理的测量值）
innovation_mahalanobis_distance = sqrt(y' / S * y);
max_acceptable_distance = 3.0;  % 3倍标准差

if innovation_mahalanobis_distance > max_acceptable_distance
    fprintf('警告: 检测到异常测量值，马氏距离 = %.2f > %.2f，降低权重处理\n', ...
            innovation_mahalanobis_distance, max_acceptable_distance);
    
    % 检查具体是哪个测量值异常
    normalized_innovations = abs(y) ./ sqrt(diag(S));
    for i = 1:length(y)
        if normalized_innovations(i) > 3.0
            fprintf('  异常测量分量 %d: 偏差 = %.2f, 标准化值 = %.2f\n', ...
                    i, y(i), normalized_innovations(i));
            
            % 自适应调整增益（降低异常分量的权重）
            adjustment_factor = 3.0 / normalized_innovations(i);
            K(:, i) = K(:, i) * adjustment_factor;
        end
    end
end

% 自适应调整卡尔曼增益，增强响应性
% 如果观测和预测差异较大，增大增益以更快响应变化
if abs(y(1)) > 10.0  % 距离差异超过10米
    K(:, 1) = K(:, 1) * 1.5;  % 增大50%
    fprintf('距离差异较大(%.2f m)，增强响应\n', y(1));
end

% 更新状态估计
kf.x = kf.x + K * y;

% 更新状态协方差（Joseph形式，数值更稳定）
I = eye(size(kf.P));
kf.P = (I - K * kf.H) * kf.P * (I - K * kf.H)' + K * kf.R * K';

% 打印更新信息
fprintf('卡尔曼更新: 预测=[%.2f, %.2f, %.2f], 测量=[%.2f, %.2f, %.2f], 残差=[%.2f, %.2f, %.2f]\n', ...
        z_pred(1), z_pred(2), z_pred(3), z(1), z(2), z(3), y(1), y(2), y(3));

end 