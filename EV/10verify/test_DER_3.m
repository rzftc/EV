%% test_DER_3.m (分布式光伏 PV 调节潜力仿真)
% 输出文件: DER3.mat

clc; clear; close all;
rng(2024);

%% 1. 初始化参数
% 假设聚合了一个小型光伏电站或光伏集群
N_PV = 1; % 聚合单元数量 (此处简化为1个大容量集群)
P_N_Cluster = 1000; % 集群额定功率 (kW), 10MW
Eta_PV = 0.95; % 逆变器及系统效率

%% 2. 时间参数定义 (与EV/AC保持一致)
dt_short = 5;     % 短时间步长 (分钟)
dt_long = 60;     % 长时间步长 (分钟)
simulation_start_hour = 6; % 仿真从早上6点开始
t_sim_hours = 24; % 仿真时长 24小时
dt = dt_short / 60; % 小时步长

% 时间轴 (小时)
time_points_absolute = simulation_start_hour : dt : (simulation_start_hour + t_sim_hours - dt);
total_steps = length(time_points_absolute);

%% 3. 创建结果存储结构
results = struct(...
    'time_points_absolute', time_points_absolute, ...
    'P_base_agg', zeros(1, total_steps), ... % 聚合基线功率 (发电为正)
    'PV_Up',      zeros(1, total_steps), ... % 上调潜力 (增加发电/减少负荷)
    'PV_Down',    zeros(1, total_steps) ...  % 下调潜力 (减少发电/增加负荷)
);

%% 4. 计算光伏基线功率 (基于 Beta 分布模拟光照)
% 技术报告公式 (3-63) 至 (3-70)
fprintf('正在计算光伏基线功率...\n');

for t = 1:total_steps
    current_h = mod(time_points_absolute(t), 24);
    
    % 简单模拟光照强度 (6:00 - 19:00 有光)
    if current_h >= 6 && current_h <= 19
        % 使用正弦波模拟理想光照基线，叠加随机波动
        % 归一化时间 x: 0 (6:00) -> 1 (19:00)
        x_norm = (current_h - 6) / (19 - 6);
        irradiance_base = sin(x_norm * pi); 
        
        % 添加随机云层遮挡波动 (Beta分布特性模拟)
        fluctuation = 0.1 * randn(); 
        irradiance = max(0, min(1, irradiance_base + fluctuation));
    else
        irradiance = 0;
    end
    
    % 计算输出功率 P = Area * eta * Irradiance -> 简化为 P_rated * Ratio
    P_gen = P_N_Cluster * irradiance * Eta_PV;
    
    results.P_base_agg(t) = P_gen;
end

%% 5. 计算调节潜力
% 调节能力定义:
% 对于电源侧: 
%   上调 (Up): 增加向电网注入功率。由于PV通常运行在MPPT，上调能力通常为0 (除非预留备用)。
%   下调 (Down): 减少向电网注入功率 (弃光)。能力为当前发电量。
% 注意: 为了与负荷侧(EV/AC)统一符号 (正值为调节能力大小):
%   Up: 增加净负荷/减少发电 (弃光) -> 这里通常称为下调(Down-regulation)对于电源，
%       但对于"与电网互济"，我们定义 Up 为 "帮助电网频率上升/弥补缺额" (增加发电)，Down 为 "帮助频率下降/消纳过剩" (减少发电)。

fprintf('正在计算光伏调节潜力...\n');

for t = 1:total_steps
    P_curr = results.P_base_agg(t);
    
    % 上调潜力 (增加出力): 
    % 假设未预留备用，MPPT模式下上调潜力为 0
    results.PV_Up(t) = 0; 
    
    % 下调潜力 (减少出力):
    % 最大可减少至 0，因此潜力大小为 P_curr
    % 符号约定: 潜力值通常为正数，表示可调节的幅度
    % 在联合仿真中，需注意方向: 减少发电 = 净负荷增加 = 消纳能力
    results.PV_Down(t) = -P_curr; % 负号表示功率减少的方向
end

%% 6. 结果保存与可视化
outputFileName = 'DER3.mat';
fprintf('正在保存结果到 %s ...\n', outputFileName);
save(outputFileName, 'results', '-v7.3');

figure('Name', 'DER3: 分布式光伏调节潜力');
plot(time_points_absolute, results.P_base_agg, 'g-', 'LineWidth', 1.5, 'DisplayName', '基线功率 (Generation)');
hold on;
plot(time_points_absolute, results.P_base_agg + results.PV_Up, 'r--', 'DisplayName', '上调边界');
plot(time_points_absolute, results.P_base_agg + results.PV_Down, 'b--', 'DisplayName', '下调边界 (弃光)');
xlabel('时间 (小时)'); ylabel('功率 (kW)');
title('分布式光伏 (PV) 时空互济潜力');
legend; grid on;