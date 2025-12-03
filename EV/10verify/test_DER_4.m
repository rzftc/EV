%% test_DER_4.m (电制氢 P2G 调节潜力仿真)
% 理论依据: 技术报告 3.2.5 电制氢模型
% 输出文件: DER4.mat

clc; clear; close all;
rng(2024);

%% 1. 初始化参数 (参考报告表 3.1)
% 假设聚合了几个典型的电制氢站
P_rated_total = 1000; % 总额定功率 20MW
P_min_maintain = 0.2 * P_rated_total; % 最小维持功率 20%
Eta_P2G = 0.75; % 电-氢转换效率
Tank_Capacity_kWh = P_rated_total * 4; % 储氢罐容量 (等效kWh)，假设能全功率运行4小时
SOC_min = 0.1;
SOC_max = 0.9;
SOC_init = 0.5;

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
    'P_base_agg', zeros(1, total_steps), ... % 基线负荷 (用电为正)
    'P2G_Up',     zeros(1, total_steps), ... % 上调潜力 (减少用电)
    'P2G_Down',   zeros(1, total_steps), ... % 下调潜力 (增加用电)
    'SOC_Tank',   zeros(1, total_steps) ...  % 储氢罐状态
);

%% 4. 生成基线功率 (工业用氢需求)
% 假设化工厂/加氢站有固定的用氢计划，折算为电力基线
% 模拟一个较为平稳但有波动的基线，白天稍高
baseline_profile = 0.6 + 0.2 * sin(2*pi*(time_points_absolute-8)/24); 
results.P_base_agg = P_rated_total * baseline_profile;

%% 5. 计算调节潜力 (考虑储罐约束)
% 技术报告公式 (3-81) 至 (3-85)

current_SOC = SOC_init;
current_E = current_SOC * Tank_Capacity_kWh;

fprintf('正在计算电制氢调节潜力...\n');

for t = 1:total_steps
    P_base = results.P_base_agg(t);
    
    % 1. 更新当前状态 (假设基线运行下的SOC变化)
    % 假设消耗氢气的速率与基线生产速率平衡 (即基线是为了满足用氢)，
    % 但为了展示灵活性，我们假设有一个独立的用氢流出 V_out
    % 简化：基线状态下 SOC 保持相对稳定 (产 = 销)
    % 这里的潜力是相对于基线的偏离能力
    
    results.SOC_Tank(t) = current_E / Tank_Capacity_kWh;
    
    % 调节时长 (假设电网指令持续 1 小时用于评估容量)
    T_regulate = 1; 
    
    % --- 上调潜力 (削减负荷 / 减少功率) ---
    % 限制 1: 功率不能低于 P_min
    P_up_power_limit = max(0, P_base - P_min_maintain);
    
    % 限制 2: 减少功率会导致产氢减少，储罐液位下降。不能低于 SOC_min
    % E_end = E_curr + (P_base - P_up)*dt*eta - V_out*dt
    % 假设基线时 dE/dt = 0 (产销平衡)，则减载导致 dE = -P_up * dt * eta
    % 允许的最大能量亏空: E_curr - E_min
    E_down_allowable = current_E - SOC_min * Tank_Capacity_kWh;
    P_up_energy_limit = E_down_allowable / (Eta_P2G * T_regulate);
    
    results.P2G_Up(t) = min(P_up_power_limit, P_up_energy_limit);
    
    % --- 下调潜力 (增加负荷 / 增加功率) ---
    % 限制 1: 功率不能高于 P_rated
    P_down_power_limit = max(0, P_rated_total - P_base);
    
    % 限制 2: 增加功率导致产氢增加，储罐液位上升。不能高于 SOC_max
    E_up_allowable = SOC_max * Tank_Capacity_kWh - current_E;
    P_down_energy_limit = E_up_allowable / (Eta_P2G * T_regulate);
    
    results.P2G_Down(t) = -min(P_down_power_limit, P_down_energy_limit); % 负号表示增加负荷
    
    % 简单的状态更新 (模拟基线下的微小波动，防止SOC一直不变)
    simulation_noise = 0.05 * P_rated_total * randn() * dt * Eta_P2G;
    current_E = max(SOC_min*Tank_Capacity_kWh, min(SOC_max*Tank_Capacity_kWh, current_E + simulation_noise));
end

%% 6. 保存与可视化
outputFileName = 'DER4.mat';
fprintf('正在保存结果到 %s ...\n', outputFileName);
save(outputFileName, 'results', '-v7.3');

figure('Name', 'DER4: 电制氢调节潜力');
subplot(2,1,1);
plot(time_points_absolute, results.P_base_agg, 'k-', 'LineWidth', 1.5, 'DisplayName', '基线功率');
hold on;
plot(time_points_absolute, results.P_base_agg - results.P2G_Up, 'r--', 'DisplayName', '最小功率边界 (上调)');
plot(time_points_absolute, results.P_base_agg - results.P2G_Down, 'b--', 'DisplayName', '最大功率边界 (下调)');
xlabel('时间 (小时)'); ylabel('功率 (kW)'); title('电制氢 (P2G) 功率调节范围'); legend; grid on;

subplot(2,1,2);
plot(time_points_absolute, results.SOC_Tank, 'm-', 'LineWidth', 1.5);
xlabel('时间 (小时)'); ylabel('储罐 SOC'); title('储氢罐状态 (基线)'); grid on;
ylim([0 1]);