%% 山东省数据中心集群能耗仿真模型 (基于2025/2030预测数据)
% 基于文档3.85-3.96公式修改
% 结合山东省数据中心负荷总规模预测进行宏观建模
% 修改时间：2025-12-08

%% 1. 清理工作区
clear; clc; close all;

%% 2. 基础物理模型参数设置 (保持单体模型物理特性不变)
% 基本参数
P_idle = 200;               % 单台服务器待机功率(W)
P_peak = 400;               % 单台服务器峰值功率(W)
mu = 1;                     % 工作负载处理速率
Delta_t = 1;                % 时间间隔1小时

% 热力学与环境参数
rho = 1.225;                % 空气密度(kg/m³)
f = 25;                     % 环境温度(℃)
kappa = 1005;               % 空气比热容(J/(kg·K))
M_2 = 0.05;                 % 热分布系数
E_red = 25;                 % 机房温度上限(℃)

% 配电参数
t_1 = 0.5;                  % 服务器数量系数
t_2 = 0.02;                 % 工作负载系数

%% 3. 规模预测目标设定 (用户输入数据)
% 目标：反推所需的等效服务器数量以匹配总负荷规模
Target_Load_2025_kW = 260 * 10^4;  % 2025年估计总规模: 260万千瓦
Target_Load_2030_kW = 900 * 10^4;  % 2030年估计总规模: 900万千瓦

%% 4. 中间参数与单机峰值能耗计算
% 计算中间变量 M1
M_1 = 1/(rho*f*kappa - M_2*rho*f*kappa) - 1/(rho*f*kappa);

% 计算送风温度 E_sup (假设按峰值功率设计冷却)
E_sup = E_red - M_1*(P_peak/1000); 

% 计算制冷效率 CoP
CoP = 0.0068*E_sup^2 + 0.0008*E_sup + 0.458;

% 计算配电系数 t_3
t_3 = t_1 / P_idle;         

% 计算综合能耗系数 nu (即理论PUE值)
% P_total = nu * P_server
nu = 1 + t_3 + 1/CoP;

% 计算"单台服务器单元"在满载情况下的总功耗 (含冷却和配电损耗)
% 单位: W
P_unit_total_peak_W = P_peak * nu; 

% 打印单机参数供参考
fprintf('---------------- 模型参数校验 ----------------\n');
fprintf('计算得出的 PUE (nu) 值: %.4f\n', nu);
fprintf('单机柜/单元满载总功耗: %.2f W\n', P_unit_total_peak_W);
fprintf('--------------------------------------------\n');

%% 5. 反推服务器规模 (S_i)
% 根据目标峰值负荷反推需要的服务器总数
% S_i = 目标总功率(W) / 单机满载总功耗(W)
S_i_2025 = (Target_Load_2025_kW * 1000) / P_unit_total_peak_W;
S_i_2030 = (Target_Load_2030_kW * 1000) / P_unit_total_peak_W;

fprintf('2025年模拟等效服务器数量: %.0f 台\n', S_i_2025);
fprintf('2030年模拟等效服务器数量: %.0f 台\n', S_i_2030);

%% 6. 24小时仿真计算函数
% 定义时间轴
T = 24; 
t_axis = 1:T;
% 定义昼夜负载波动曲线 (0.6 ~ 1.0 波动率)
% 模拟全省数据中心并非时刻满载，存在夜间低谷
load_profile = 0.8 + 0.2*sin(2*pi*(t_axis-10)/24); 

% 封装计算过程
[P_total_2025, P_comp_2025] = calc_load_profile(S_i_2025, P_idle, P_peak, mu, Delta_t, CoP, t_3, load_profile, T);
[P_total_2030, P_comp_2030] = calc_load_profile(S_i_2030, P_idle, P_peak, mu, Delta_t, CoP, t_3, load_profile, T);

% 单位转换 W -> 万kW
P_total_2025_wkW = P_total_2025 / 10^7;
P_total_2030_wkW = P_total_2030 / 10^7;

%% 7. 绘图：总负荷对比 (图1)
figureConfig = {'Color','white', 'Position', [100 100 1200 900]};
fig1 = figure(figureConfig{:});
set(gca, 'FontName', 'Microsoft YaHei', 'FontSize', 14);

plot(t_axis, P_total_2030_wkW, 'r-o', 'LineWidth', 2.5, 'MarkerSize', 8, 'DisplayName', '2030年预测 (900万kW级)');
hold on;
plot(t_axis, P_total_2025_wkW, 'b-s', 'LineWidth', 2.5, 'MarkerSize', 8, 'DisplayName', '2025年预测 (260万kW级)');

grid on;
xlabel('时间 (小时)', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('总负荷 (万千瓦)', 'FontSize', 16, 'FontWeight', 'bold');
title('山东省数据中心集群时序负荷预测 (2025 vs 2030)', 'FontSize', 18, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 14);
xlim([1 24]);
xticks(1:2:24);

% 添加峰值标注
[max_2030, idx_2030] = max(P_total_2030_wkW);
text(idx_2030, max_2030, sprintf('  峰值: %.1f 万kW', max_2030), 'FontSize', 12, 'Color', 'red', 'FontWeight', 'bold');

[max_2025, idx_2025] = max(P_total_2025_wkW);
text(idx_2025, max_2025, sprintf('  峰值: %.1f 万kW', max_2025), 'FontSize', 12, 'Color', 'blue', 'FontWeight', 'bold');

print(fig1, 'SD_DataCenter_Projection_Comparison.png', '-dpng', '-r300');

%% 8. 绘图：2030年能耗分解 (图2)
% 展示900万千瓦规模下的内部构成
fig2 = figure(figureConfig{:});
set(gca, 'FontName', 'Microsoft YaHei', 'FontSize', 14);

% 准备绘图数据 (单位: 万kW)
p_server = P_comp_2030.servers / 10^7;
p_cool   = P_comp_2030.cool / 10^7;
p_net    = P_comp_2030.net / 10^7;

area_data = [p_server; p_cool; p_net]';
h = area(t_axis, area_data);

% 样式调整
h(1).FaceColor = [0.2 0.6 0.8]; % 服务器
h(2).FaceColor = [0.4 0.8 0.4]; % 冷却
h(3).FaceColor = [0.9 0.7 0.3]; % 配电

grid on;
xlabel('时间 (小时)', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('功率 (万千瓦)', 'FontSize', 16, 'FontWeight', 'bold');
title('2030年山东省数据中心负荷成分分解 (堆积图)', 'FontSize', 18, 'FontWeight', 'bold');
legend({'IT设备能耗', '冷却系统能耗', '配电系统损耗'}, 'Location', 'northwest', 'FontSize', 14);
xlim([1 24]);
xticks(1:2:24);
ylim([0 1000]); % 略高于900以便留白

print(fig2, 'SD_DataCenter_2030_Breakdown.png', '-dpng', '-r300');

%% 辅助函数：负荷计算内核
function [P_total, P_comp] = calc_load_profile(S_i, P_idle, P_peak, mu, Delta_t, CoP, t_3, load_ratio, T)
    P_servers = zeros(1,T);
    P_cool = zeros(1,T);
    P_network = zeros(1,T);
    P_total = zeros(1,T);
    
    for k = 1:T
        % 当前时刻的单机负载 (在P_idle和P_peak之间波动)
        % load_ratio(k) 是 0~1 之间的负载率，实际上公式3.85用的是任务量lambda
        % 这里简化为：实际功率 = 待机 + (峰值-待机)*负载率
        P_server_unit = P_idle + (P_peak - P_idle) * load_ratio(k);
        
        % 1. 服务器集群总能耗
        P_servers(k) = S_i * P_server_unit;
        
        % 2. 冷却系统能耗 (公式3.92: P_cool = P_server / CoP)
        P_cool(k) = P_servers(k) / CoP;
        
        % 3. 配电系统能耗 (公式3.93简化: P_net = t_3 * P_server)
        P_network(k) = t_3 * P_servers(k);
        
        % 4. 总能耗
        P_total(k) = P_servers(k) + P_cool(k) + P_network(k);
    end
    
    % 返回分项数据结构
    P_comp.servers = P_servers;
    P_comp.cool = P_cool;
    P_comp.net = P_network;
end