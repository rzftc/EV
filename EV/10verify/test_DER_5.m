%% test_DER_5.m (数据中心 IDC 调节潜力仿真)
% 理论依据: 技术报告 3.2.6 数据中心模型
% 输出文件: DER5.mat

clc; clear; close all;
rng(2024);

%% 1. 初始化参数
% 假设一个中型数据中心
N_Server = 1000; % 服务器数量
P_server_idle = 0.15; % 待机功率 kW
P_server_peak = 0.40; % 峰值功率 kW
PUE = 1.4; % 电能利用效率 (Power Usage Effectiveness)

%% 2. 时间参数
dt_short = 5;     
simulation_start_hour = 6; 
t_sim_hours = 24; 
dt = dt_short / 60; 
time_points_absolute = simulation_start_hour : dt : (simulation_start_hour + t_sim_hours - dt);
total_steps = length(time_points_absolute);

%% 3. 创建结果存储
results = struct(...
    'time_points_absolute', time_points_absolute, ...
    'P_base_agg', zeros(1, total_steps), ... 
    'IDC_Up',     zeros(1, total_steps), ... 
    'IDC_Down',   zeros(1, total_steps) ...  
);

%% 4. 生成基线功率 (IT 负载 + 冷却)
% 技术报告公式 (3-86) 至 (3-97)
fprintf('正在计算数据中心基线功率...\n');

% 模拟 IT 负载率 (Workload Rate, 0-1)
% 通常白天繁忙，夜间有批处理任务
workload_profile = 0.5 + 0.3 * sin(2*pi*(time_points_absolute-9)/24); 
% 添加随机性
workload_profile = max(0.1, min(1.0, workload_profile + 0.05*randn(1, total_steps)));

% 计算 P_IT (服务器功率)
% P_server = N * (P_idle + u * (P_peak - P_idle))
P_IT = N_Server * (P_server_idle + workload_profile .* (P_server_peak - P_server_idle));

% 计算总功率 P_total = P_IT * PUE
results.P_base_agg = P_IT * PUE;

%% 5. 计算调节潜力
% 调节来源:
% 1. 算力负荷转移 (Delay-tolerant jobs): 假设 20% 的负载是可延迟的批处理任务
% 2. 冷却系统调节 (短期): 利用热惯性，可短期减少 10% 冷却功率或预冷增加 10%

ratio_flexible_IT = 0.20; % 20% IT 负载可调节
ratio_flexible_Cooling = 0.15; % 15% 冷却负载可调节 (P_cooling = P_total - P_IT)

fprintf('正在计算数据中心调节潜力...\n');

for t = 1:total_steps
    P_base = results.P_base_agg(t);
    P_it_curr = P_IT(t);
    P_cooling_curr = P_base - P_it_curr;
    
    % 上调潜力 (削减负荷):
    % 1. 延迟批处理任务: 减少 P_IT
    delta_P_IT_up = P_it_curr * ratio_flexible_IT; 
    % 2. 提高空调温度设定: 减少 P_Cooling
    delta_P_Cool_up = P_cooling_curr * ratio_flexible_Cooling;
    
    % 总上调潜力 (考虑 PUE 耦合，简化为各部分之和)
    results.IDC_Up(t) = delta_P_IT_up * PUE; % 假设 IT 减少带动冷却减少，近似用 PUE 放大
    
    % 下调潜力 (增加负荷):
    % 1. 提前执行任务 / 增加算力利用: 增加 P_IT (受限于服务器峰值)
    P_IT_max = N_Server * P_server_peak;
    delta_P_IT_down = min(P_IT_max - P_it_curr, P_it_curr * ratio_flexible_IT);
    % 2. 预冷: 增加 P_Cooling
    delta_P_Cool_down = P_cooling_curr * ratio_flexible_Cooling;
    
    results.IDC_Down(t) = -(delta_P_IT_down * PUE); % 负号
end

%% 6. 保存与可视化
outputFileName = 'DER5.mat';
fprintf('正在保存结果到 %s ...\n', outputFileName);
save(outputFileName, 'results', '-v7.3');

figure('Name', 'DER5: 数据中心调节潜力');
plot(time_points_absolute, results.P_base_agg, 'k-', 'LineWidth', 1.5, 'DisplayName', '基线总功率');
hold on;
fill([time_points_absolute fliplr(time_points_absolute)], ...
     [results.P_base_agg-results.IDC_Up fliplr(results.P_base_agg-results.IDC_Down)], ...
     'c', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', '调节可行域');
xlabel('时间 (小时)'); ylabel('功率 (kW)');
title('数据中心 (IDC) 算力与温控联合调节潜力');
legend; grid on;