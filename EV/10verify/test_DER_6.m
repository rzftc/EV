%% test_DER_6.m (工业大负荷调节潜力仿真)
% 理论依据: 技术报告 3.2.7 工业大负荷模型
% 输出文件: DER6.mat

clc; clear; close all;
rng(2024);

%% 1. 初始化参数
% 模拟一个钢铁厂或大型制造厂
% 负荷组成:
% A. 关键生产线 (不可中断，但可平移班次): 功率 5000 kW
% B. 辅助设备 (风机/水泵，可连续调节): 功率 1000 kW
P_Production_Rated = 5000;
P_Aux_Rated = 1000;

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
    'Ind_Up',     zeros(1, total_steps), ... 
    'Ind_Down',   zeros(1, total_steps) ...  
);

%% 4. 生成基线功率 (排班表)
% 技术报告公式 (3-98) 至 (3-101)
fprintf('正在计算工业大负荷基线功率...\n');

P_base_prod = zeros(1, total_steps);
P_base_aux = zeros(1, total_steps);

for t = 1:total_steps
    current_h = mod(time_points_absolute(t), 24);
    
    % 生产排班: 8:00-12:00, 14:00-18:00, 20:00-04:00 (三班倒，中间有休息)
    is_working = (current_h >= 8 && current_h < 12) || ...
                 (current_h >= 14 && current_h < 18) || ...
                 (current_h >= 20 || current_h < 4);
    
    if is_working
        % 生产负荷 (加一点随机波动模拟冲击特性)
        P_base_prod(t) = P_Production_Rated * (0.9 + 0.2*rand());
    else
        P_base_prod(t) = 0; % 停机维护/休息
    end
    
    % 辅助负荷: 持续运行，白天稍高
    P_base_aux(t) = P_Aux_Rated * (0.8 + 0.2 * sin(pi*(current_h-6)/12)); 
    P_base_aux(t) = max(0, P_base_aux(t));
end

results.P_base_agg = P_base_prod + P_base_aux;

%% 5. 计算调节潜力
% 调节来源:
% 1. 生产负荷: 
%    - 上调 (中断/避峰): 仅在生产时段有效。假设可中断/降功率 50% (减产运行)。
%    - 下调 (填谷/加班): 在非生产时段可启动生产 (假设有备用产能)。
% 2. 辅助负荷: 
%    - 连续调节 +/- 20%

fprintf('正在计算工业负荷调节潜力...\n');

for t = 1:total_steps
    p_prod = P_base_prod(t);
    p_aux = P_base_aux(t);
    
    % --- 上调潜力 (削减负荷) ---
    % 辅助设备可减 20%
    d_aux_up = p_aux * 0.2;
    % 生产设备: 如果在运行，假设可降低 30% 负荷 (降速/单线运行)
    if p_prod > 100 % 正在运行
        d_prod_up = p_prod * 0.3;
    else
        d_prod_up = 0;
    end
    results.Ind_Up(t) = d_aux_up + d_prod_up;
    
    % --- 下调潜力 (增加负荷) ---
    % 辅助设备可增 20%
    d_aux_down = p_aux * 0.2;
    % 生产设备: 
    % 如果未运行(休息时段)，可启动生产 (填谷) -> 潜力巨大 (0 -> Rated)
    % 如果在运行，假设已接近满载，增加空间较小 (e.g. 10%)
    if p_prod > 100 
        d_prod_down = p_prod * 0.1;
    else
        % 休息时段，具备全额启动能力 (假设为了消纳新能源可临时开工)
        d_prod_down = P_Production_Rated;
    end
    results.Ind_Down(t) = -(d_aux_down + d_prod_down);
end

%% 6. 保存与可视化
outputFileName = 'DER6.mat';
fprintf('正在保存结果到 %s ...\n', outputFileName);
save(outputFileName, 'results', '-v7.3');

figure('Name', 'DER6: 工业大负荷调节潜力');
subplot(2,1,1);
area(time_points_absolute, results.P_base_agg, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'k');
hold on;
plot(time_points_absolute, results.P_base_agg - results.Ind_Up, 'r', 'LineWidth', 2, 'DisplayName', '最大削峰 (Up)');
ylabel('功率 (kW)'); title('工业负荷基线与削峰能力'); legend; grid on;

subplot(2,1,2);
plot(time_points_absolute, results.Ind_Down, 'b', 'LineWidth', 2);
ylabel('下调潜力 (kW)'); xlabel('时间 (小时)'); 
title('工业负荷填谷能力 (负值代表增加负荷)'); grid on;