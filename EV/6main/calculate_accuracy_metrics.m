%% calculate_accuracy_metrics.m
% 用于衡量“单体求和功率”与“聚合模型功率”误差度的主程序
clc; clear; close all;

%% 1. 加载数据
resultsFile = 'main_potential_5min_1000_8am_bound.mat'; % 确保文件名与您保存的一致
if ~exist(resultsFile, 'file')
    error('未找到结果文件: %s，请先运行主仿真程序。', resultsFile);
end
fprintf('正在加载结果文件: %s ...\n', resultsFile);
data = load(resultsFile);
results = data.results;

%% 2. 提取数据
% P_true: 实际统计值（基准），即 results.P_agg
% P_pred: 聚合模型计算值，即 results.EV_Power
if ~isfield(results, 'EV_Power')
    error('结果结构体中缺少 EV_Power 字段，请检查是否使用了最新版的主程序。');
end

P_true = results.P_agg;      % 真实值 (Benchmark)
P_pred = results.EV_Power;   % 预测值 (Model)

% 确保数据为列向量以便统一计算
P_true = P_true(:);
P_pred = P_pred(:);

%% 3. 计算误差指标

% (1) 绝对误差向量
abs_error = abs(P_true - P_pred);

% (2) 均方根误差 (RMSE - Root Mean Square Error)
% 反映误差的离散程度，值越小越好
rmse_val = sqrt(mean((P_true - P_pred).^2));

% (3) 平均绝对误差 (MAE - Mean Absolute Error)
% 反映误差的平均水平
mae_val = mean(abs_error);

% (4) 最大绝对误差 (Max Error)
% 反映最坏情况下的偏差
[max_error_val, max_error_idx] = max(abs_error);

% (5) 平均绝对百分比误差 (MAPE - Mean Absolute Percentage Error)
% 仅在 P_true 不为 0 (或大于微小阈值) 的时刻计算，避免除以零
threshold = 1e-3; % 1瓦特以下的功率视为0
valid_mask = abs(P_true) > threshold;

if any(valid_mask)
    mape_val = mean( abs((P_true(valid_mask) - P_pred(valid_mask)) ./ P_true(valid_mask)) ) * 100;
else
    mape_val = 0;
end

% (6) 决定系数 (R-Squared)
% 衡量拟合优度，越接近 1 越好
SS_res = sum((P_true - P_pred).^2);       % 残差平方和
SS_tot = sum((P_true - mean(P_true)).^2); % 总平方和
r_squared = 1 - (SS_res / SS_tot);

%% 4. 输出结果报告
fprintf('\n========================================\n');
fprintf('       聚合模型精度评估报告       \n');
fprintf('========================================\n');
fprintf('1. 均方根误差 (RMSE)    : %.4f kW\n', rmse_val);
fprintf('2. 平均绝对误差 (MAE)   : %.4f kW\n', mae_val);
fprintf('3. 最大绝对误差 (MaxErr): %.4f kW (发生在第 %d 步)\n', max_error_val, max_error_idx);
fprintf('4. 平均百分比误差 (MAPE): %.2f %%\n', mape_val);
fprintf('5. 拟合优度 (R-Squared) : %.4f\n', r_squared);
fprintf('----------------------------------------\n');

% 简单的判定
if mape_val < 5 && r_squared > 0.95
    fprintf('>> 评估结论: 模型精度 [优秀]\n');
elseif mape_val < 10
    fprintf('>> 评估结论: 模型精度 [良好]\n');
else
    fprintf('>> 评估结论: 模型精度 [一般]，建议检查分组逻辑\n');
end
fprintf('========================================\n');

%% 5. (可选) 绘制误差分布直方图
figure('Name', '误差分布分析', 'Position', [300 300 600 400], 'NumberTitle', 'off');
histogram(P_true - P_pred, 50, 'FaceColor', [0.2 0.6 0.8]);
xlabel('误差值 (P_{agg} - EV_{Power}) [kW]');
ylabel('频数');
title(['误差分布直方图 (RMSE = ' num2str(rmse_val, '%.2f') ' kW)']);
grid on;