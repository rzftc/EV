%% 山东省分布式储能(ESS)时空互济潜力预测主程序 (修正版)
% 基于技术报告第3章分布式储能聚合模型 (公式3.52 - 3.61)
% 数据来源: 用户提供山东省2024/2030规划数据
% 修正内容: 修复了 calc_ess_potential 函数中字段名 E_rated_MW 的拼写错误
% 时间: 2025-12-08

clear; clc; close all;

%% 1. 基础参数设置 (调研数据)
% 典型单体ESS参数
ESS_Unit.P_rated_kW = 186;       % 单机额定功率 (kW)
ESS_Unit.E_rated_kWh = 372;      % 单机额定容量 (kWh)
ESS_Unit.P_rated_MW = ESS_Unit.P_rated_kW / 1000; % 转换为MW
ESS_Unit.E_rated_MWh = ESS_Unit.E_rated_kWh / 1000; % 转换为MWh (注意这里定义的是 MWh)
ESS_Unit.eta_ch = 0.95;          % 充电效率 (假设)
ESS_Unit.eta_dis = 0.95;         % 放电效率 (假设)
ESS_Unit.SOC_min = 0.1;          % 最小SOC约束
ESS_Unit.SOC_max = 0.9;          % 最大SOC约束

% 规模场景数据 (MW)
Scale_2024_MW = 28;              % 2024年分布式储能总规模
Scale_2030_MW = 182;             % 2030年分布式储能总规模预计

% 仿真设置
T = 24;                          % 仿真时长 (小时)
dt = 1;                          % 时间步长 (小时)
time_axis = 1:T;

%% 2. 规模反推与初始化
% 计算等效聚合台数
Num_Units_2024 = round(Scale_2024_MW / ESS_Unit.P_rated_MW);
Num_Units_2030 = round(Scale_2030_MW / ESS_Unit.P_rated_MW);

fprintf('---------- 分布式储能聚合规模概览 ----------\n');
fprintf('单体参数: P=%.3f MW, E=%.3f MWh (2小时系统)\n', ESS_Unit.P_rated_MW, ESS_Unit.E_rated_MWh);
fprintf('[2024年] 总规模: %d MW, 等效聚合台数: %d 台\n', Scale_2024_MW, Num_Units_2024);
fprintf('[2030年] 总规模: %d MW, 等效聚合台数: %d 台\n', Scale_2030_MW, Num_Units_2030);
fprintf('--------------------------------------------\n');

%% 3. 聚合模型仿真 (计算调节潜力)
% 假设一个基础运行场景：储能系统在一天中维持一定SOC水平(如50%)备用，
% 或者随时间发生自适应变化。此处计算其在任意时刻的"最大可调用潜力"。

% 场景1: 2024年潜力预测
[P_up_2024, P_down_2024] = calc_ess_potential(Num_Units_2024, ESS_Unit, T, dt);

% 场景2: 2030年潜力预测
[P_up_2030, P_down_2030] = calc_ess_potential(Num_Units_2030, ESS_Unit, T, dt);

%% 4. 结果可视化
figure('Color', 'white', 'Position', [100, 100, 1200, 600]);

% 子图1: 2024年调节潜力
subplot(1, 2, 1);
fill([time_axis, fliplr(time_axis)], [P_up_2024, fliplr(zeros(size(P_up_2024)))], [0.6, 0.8, 1], 'EdgeColor', 'none', 'FaceAlpha', 0.5); hold on;
fill([time_axis, fliplr(time_axis)], [P_down_2024, fliplr(zeros(size(P_down_2024)))], [1, 0.6, 0.6], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
plot(time_axis, P_up_2024, 'b-', 'LineWidth', 2);
plot(time_axis, P_down_2024, 'r-', 'LineWidth', 2);
yline(0, 'k-', 'LineWidth', 1);
title('2024年山东省分布式储能调节潜力预测 (总规模28MW)');
xlabel('时间 (小时)');
ylabel('调节功率 (MW)');
legend('上调节潜力(放电)', '下调节潜力(充电)', 'Location', 'best');
grid on;
% ylim根据数据范围自动调整或手动设置
% ylim([-40 40]); 

% 子图2: 2030年调节潜力
subplot(1, 2, 2);
fill([time_axis, fliplr(time_axis)], [P_up_2030, fliplr(zeros(size(P_up_2030)))], [0.6, 0.8, 1], 'EdgeColor', 'none', 'FaceAlpha', 0.5); hold on;
fill([time_axis, fliplr(time_axis)], [P_down_2030, fliplr(zeros(size(P_down_2030)))], [1, 0.6, 0.6], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
plot(time_axis, P_up_2030, 'b-', 'LineWidth', 2);
plot(time_axis, P_down_2030, 'r-', 'LineWidth', 2);
yline(0, 'k-', 'LineWidth', 1);
title('2030年山东省分布式储能调节潜力预测 (总规模182MW)');
xlabel('时间 (小时)');
ylabel('调节功率 (MW)');
legend('上调节潜力(放电)', '下调节潜力(充电)', 'Location', 'best');
grid on;
% ylim([-250 250]); 

%% 5. 核心计算函数
function [P_up_max, P_down_max] = calc_ess_potential(N, Unit, T, dt)
    % 输入: 
    % N: 聚合单体数量
    % Unit: 单体参数结构体
    % T: 时间长度
    % dt: 时间步长
    
    % 初始化输出
    P_up_max = zeros(1, T);   % 最大上调节(放电, >0)
    P_down_max = zeros(1, T); % 最大下调节(充电, <0)
    
    % 模拟初始SOC分布
    % 假设聚合体内部SOC服从正态分布，均值0.5 (备用状态)
    soc_mean = 0.5;
    soc_std = 0.1;
    SOC_dist = soc_mean + soc_std * randn(N, 1);
    
    % 修正越界SOC
    SOC_dist(SOC_dist > Unit.SOC_max) = Unit.SOC_max;
    SOC_dist(SOC_dist < Unit.SOC_min) = Unit.SOC_min;
    
    % 计算每一时刻的瞬时潜力
    for t = 1:T
        % 模拟SOC随时间的自然波动
        current_SOCs = SOC_dist + 0.05 * sin(2*pi*t/24); 
        current_SOCs(current_SOCs > Unit.SOC_max) = Unit.SOC_max;
        current_SOCs(current_SOCs < Unit.SOC_min) = Unit.SOC_min;
        
        % 1. 计算上调节能力 (Discharge, 减少电网负荷/反向送电)
        % 修正点：使用 E_rated_MWh
        E_current = current_SOCs .* Unit.E_rated_MWh; 
        E_min_val = Unit.SOC_min * Unit.E_rated_MWh;
        
        p_dis_available = (E_current - E_min_val) ./ dt .* Unit.eta_dis;
        p_up_units = min(Unit.P_rated_MW, p_dis_available);
        P_up_max(t) = sum(p_up_units);
        
        % 2. 计算下调节能力 (Charge, 增加负荷)
        % 修正点：使用 E_rated_MWh
        E_max_val = Unit.SOC_max * Unit.E_rated_MWh;
        
        p_ch_available = (E_max_val - E_current) ./ dt ./ Unit.eta_ch;
        p_down_units = min(Unit.P_rated_MW, p_ch_available);
        P_down_max(t) = -sum(p_down_units); % 约定下调节为负值
    end
end