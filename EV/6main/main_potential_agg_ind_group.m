%% main_potential_agg_ind_8am.m
% (新增功能: 同时计算并对比“单体求和”潜力 - 包含分组聚合改进)
% (本次修改: 增加聚合模型计算的实时功率 results.EV_Power，用于验证)
% (本次修改: 增加 results.P_agg_ptcp，保存参与聚合EV的运行功率)
% (本次修改: 特定组(PN=6.6, C=26.9, 参与聚合)的期望SOC与Lambda对比，保存高DPI图)
% (本次修改: 强制同质化参数，并修正 S_agg_prime 的计算逻辑为直接平均)

clc; clear; close all;
rng(2024);

%% 初始化参数
excelFile = 'evtest_8am.xlsx';
if ~exist(excelFile, 'file')
    generateEVParameters_real(excelFile, 1000, 1.0);
    fprintf('已生成参数模板: %s\n', excelFile);
end
[EVs, t_sim, ~, ~, P_tar] = initializeFromExcel_8am(excelFile);
fprintf('成功加载%d辆EV数据\n', length(EVs));
    
%% 时间参数定义
dt_short = 5;     % 短时间步长 (分钟)
dt_long = 60;       % 长时间步长 (分钟)
simulation_start_hour = 8;
simulation_end_hour   = 32;
dt = dt_short / 60;       % 短步长 (小时)
dt_minutes = dt_short;    % 短步长 (分钟)
dt_long_minutes = dt_long;% 长步长 (分钟)
t_adj = 5 / 60;     % 调节时长 (小时)
% P_tar=zeros(1,24);
% 时间轴 (小时)
time_points_absolute = simulation_start_hour : dt : (simulation_start_hour + t_sim/60 - dt);
num_time_points = length(time_points_absolute);
fprintf('仿真时间范围: %.2f 小时 到 %.2f 小时, 步长: %.3f 小时 (%d 点)\n', ...
    time_points_absolute(1), time_points_absolute(end), dt, num_time_points);

%% 索引预处理
num_long_steps = t_sim / dt_long;
num_short_per_long = dt_long / dt_short;
total_steps = num_long_steps * num_short_per_long;
assert(mod(num_long_steps, 24*(60/dt_long)) == 0, '长时间步数必须是24小时对应步数的整数倍');

%% 创建结果存储结构
num_evs = length(EVs);
results = struct(...
    'P_agg',         zeros(1, total_steps), ...
    'P_agg_ptcp',    zeros(1, total_steps), ... % [新增] 参与聚合的EV运行功率
    'P_base',        zeros(1, total_steps), ...
    'S_agg',         zeros(1, total_steps), ...
    'lambda',        zeros(1, total_steps), ...
    'EV_S_original', zeros(num_evs, total_steps), ...
    'EV_S_mod',      zeros(num_evs, total_steps), ...
    'm3',            zeros(num_evs, 1), ...
    'P_tar',         zeros(1, total_steps), ...
    'P_cu',          zeros(1, total_steps), ...
    'EV_Up',         zeros(1, total_steps), ... % 聚合模型上调潜力
    'EV_Down',       zeros(1, total_steps), ... % 聚合模型下调潜力
    'EV_Power',      zeros(1, total_steps), ... % [新增] 聚合模型实时功率
    'EV_Up_Individual_Sum',   zeros(1, total_steps), ... % 单体求和上调潜力
    'EV_Down_Individual_Sum', zeros(1, total_steps),  ... % 单体求和下调潜力
    ...
    'SOC_EV',             zeros(num_evs, total_steps), ... 
    'EV_Up_Individual',   zeros(num_evs, total_steps), ... 
    'EV_Down_Individual', zeros(num_evs, total_steps), ...
    ...
    'EV_E_actual',        zeros(num_evs, total_steps), ... % [新增] 实际电量轨迹
    'EV_E_baseline',      zeros(num_evs, total_steps), ... % [新增] 基线电量轨迹
    'P_base_agg',         zeros(1, total_steps), ...       % [新增] 聚合基线功率
    'EV_t_in',            zeros(num_evs, 1), ...           % [新增] 入网时间分布
    'EV_t_dep',           zeros(num_evs, 1) ...            % [新增] 离网时间分布
    );

%% 初始化EV参数 (向量化)
p_real_vec = [EVs.p_real]';
P_h_max_vec = [EVs.P_h_max]';
P_0_vec = [EVs.P_0]';
P_l_min_vec = [EVs.P_l_min]';
Delta_E_h_max_vec = [EVs.Delta_E_h_max]';
Delta_E_q_max_vec = [EVs.Delta_E_q_max]';
E_tar_set_vec = [EVs.E_tar_set]';
E_ini_vec = [EVs.E_ini]';
P_N_vec = [EVs.P_N]';
C_vec = [EVs.C]';

E_tar_original_vec = E_tar_set_vec;
SOC_original_vec = zeros(num_evs, 1);
E_reg_min_vec = E_tar_original_vec;
E_reg_max_vec = E_tar_original_vec;

p_min = 15; p_max = 50;
p_min_prime = 10; p_max_prime = 40;
base_Price_vec = 30 * ones(num_evs, 1);
p_incentive_vec = [EVs.p_incentive]';
participation_probabilities_vec = calculateParticipation(p_incentive_vec, base_Price_vec);
ptcp_vec = (rand(num_evs, 1) < participation_probabilities_vec);

E_tar_max_flex_vec = 0.2 * C_vec;
[deltaE_up_vec, deltaE_down_vec] = incentiveTempEV_updown(p_incentive_vec, p_min, p_max, p_min_prime, p_max_prime, E_tar_max_flex_vec);

participating_indices = find(ptcp_vec);
if ~isempty(participating_indices)
    E_reg_min_vec(participating_indices) = E_tar_original_vec(participating_indices) - deltaE_down_vec(participating_indices);
    E_reg_min_vec(participating_indices) = max(E_reg_min_vec(participating_indices), E_ini_vec(participating_indices));

    E_reg_max_vec(participating_indices) = E_tar_original_vec(participating_indices) + deltaE_up_vec(participating_indices);
    E_reg_max_vec(participating_indices) = min(E_reg_max_vec(participating_indices), C_vec(participating_indices));
end

% 写回EVs结构体
E_tar_orig_cell = num2cell(E_tar_original_vec);
SOC_orig_cell = num2cell(SOC_original_vec);
ptcp_cell = num2cell(ptcp_vec);
E_reg_min_cell = num2cell(E_reg_min_vec);
E_reg_max_cell = num2cell(E_reg_max_vec);
[EVs.E_tar_original] = E_tar_orig_cell{:};
[EVs.SOC_original] = SOC_orig_cell{:};
[EVs.ptcp] = ptcp_cell{:};
[EVs.E_reg_min] = E_reg_min_cell{:};
[EVs.E_reg_max] = E_reg_max_cell{:};

% 初始化 E_tar, tau_rem, E_exp, E_actual, E_current, P_current
E_tar_vec = max(E_tar_set_vec, E_ini_vec);
t_ch_vec = zeros(num_evs, 1);
valid_PN_mask = P_N_vec > 0;
if any(valid_PN_mask)
    t_ch_vec(valid_PN_mask) = 60 .* (E_tar_vec(valid_PN_mask) - E_ini_vec(valid_PN_mask)) ./ P_N_vec(valid_PN_mask);
end
t_ch_vec = max(t_ch_vec, 0);

E_tar_cell = num2cell(E_tar_vec);
[EVs.E_tar] = E_tar_cell{:};
tau_rem_cell = num2cell(t_ch_vec);
[EVs.tau_rem] = tau_rem_cell{:};
E_ini_cell = num2cell(E_ini_vec);
[EVs.E_exp] = E_ini_cell{:};
[EVs.E_actual] = E_ini_cell{:};
[EVs.E_current] = E_ini_cell{:};
P_current_init_cell = num2cell(zeros(num_evs, 1));
[EVs.P_current] = P_current_init_cell{:};

%% [新增] 定义特定观察组 (PN=6.6, C=26.9, 且参与聚合)
target_PN = 6.6;
target_C = 26.86; % 注意excel中可能存在的微小精度差异，这里保持一致
% 筛选: 物理参数匹配 AND 实际参与聚合(ptcp=1)
group_indices = find(abs([EVs.P_N] - target_PN) < 1e-3 & ...
                     abs([EVs.C] - target_C) < 1e-3 & ...
                     [EVs.ptcp] == 1);
fprintf('特定观察组 (P_N=%.1f, C=%.1f, 参与聚合) 共有 %d 辆车\n', target_PN, target_C, length(group_indices));

% [本次新增修改] 强制将该组内所有车辆的参数设置为完全一致
if ~isempty(group_indices)
    fprintf('正在将特定组 (ID: %d ...) 的 %d 辆车参数(时间/电量)设置为完全一致...\n', group_indices(1), length(group_indices));
    ref_idx = group_indices(1);
    ref_EV = EVs(ref_idx);
    
    for k = 1:length(group_indices)
        idx = group_indices(k);
        % 1. 强制统一时间参数
        EVs(idx).t_in = ref_EV.t_in;
        EVs(idx).t_dep = ref_EV.t_dep;
        
        % 2. 强制统一电量参数
        EVs(idx).E_ini = ref_EV.E_ini;
        EVs(idx).E_tar = ref_EV.E_tar;
        EVs(idx).E_tar_original = ref_EV.E_tar_original;
        EVs(idx).SOC_original = ref_EV.SOC_original;
        
        % 3. 强制统一调节边界
        EVs(idx).E_reg_min = ref_EV.E_reg_min;
        EVs(idx).E_reg_max = ref_EV.E_reg_max;
        
        % 4. 重置状态变量 (确保E_exp等从新的起点开始)
        EVs(idx).E_exp = EVs(idx).E_ini;
        EVs(idx).E_actual = EVs(idx).E_ini;
        EVs(idx).E_current = EVs(idx).E_ini;
        EVs(idx).eta = ref_EV.eta; % <--- 关键！
        EVs(idx).r   = ref_EV.r;   % <--- 关键！
        % 5. 重新计算充电视窗
        t_ch = 60 * (EVs(idx).E_tar - EVs(idx).E_ini) / EVs(idx).P_N;
        EVs(idx).tau_rem = max(t_ch, 0);
    end
end

% [新增] 在 results 结构体中增加存储空间
results.S_agg_group_raw = zeros(1, total_steps);   % 存储该组的原始聚合SOC (S)
results.S_agg_group_prime = zeros(1, total_steps); % 存储该组的变换后聚合SOC (S') -> 用于和 lambda 比较

%% 初始化聚合SOC 和 基线功率
S_agg_current = mean([EVs.SOC_original]);

fprintf('正在预计算所有EV的基线功率序列...\n');
EVs_for_baseline = EVs;
parfor i = 1:num_evs
    t_dep_h = EVs_for_baseline(i).t_dep / 60;
    t_in_h = EVs_for_baseline(i).t_in / 60;

    EVs_for_baseline(i).P_base_sequence = EVbaseP_ChargeUntilFull(...
        EVs_for_baseline(i).C, EVs_for_baseline(i).eta,...
        EVs_for_baseline(i).E_tar_original, EVs_for_baseline(i).E_ini,...
        t_dep_h, t_in_h, dt, ...
        EVs_for_baseline(i).r, EVs_for_baseline(i).P_N, ...
        EVs_for_baseline(i).SOC_original, num_time_points, time_points_absolute);
end
for i = 1:num_evs
    EVs(i).P_base_sequence = EVs_for_baseline(i).P_base_sequence;
end
clear EVs_for_baseline;
fprintf('基线功率序列计算完成。\n');

% [新增] 计算所有EV的基线能量轨迹 和 聚合基线功率
fprintf('正在计算基线能量轨迹与聚合基线功率...\n');
temp_P_base_agg = zeros(1, total_steps);
for i = 1:num_evs
    % 1. 累加聚合基线功率
    % 确保 P_base_sequence 是行向量 (1 x total_steps)
    p_base_seq = EVs(i).P_base_sequence;
    if size(p_base_seq, 1) > 1
        p_base_seq = p_base_seq'; 
    end
    temp_P_base_agg = temp_P_base_agg + p_base_seq;

    % 2. 计算单体基线能量轨迹
    % 积分: E(t) = E_ini + cumsum(P_base * eta * dt)
    delta_E = p_base_seq * EVs(i).eta * (dt_short / 60);
    results.EV_E_baseline(i, :) = EVs(i).E_ini + cumsum(delta_E);
end
results.P_base_agg = temp_P_base_agg; % 保存聚合基线功率
results.EV_t_in = [EVs.t_in]';       % 保存入网时间
results.EV_t_dep = [EVs.t_dep]';     % 保存离网时间
fprintf('基线数据计算完成。\n');


%% 外层循环（长时间步长）
for long_idx = 1:num_long_steps
    t_long_start_minute = (long_idx - 1) * dt_long_minutes;
    t_current_long_abs_minute = t_long_start_minute + simulation_start_hour * 60;
    parfor i = 1:num_evs
        EVs(i) = updateLockState(EVs(i), t_current_long_abs_minute); 
    end
    %% 长时间步处理
    [lambda_star] = aggregateEVs(EVs, P_tar(long_idx));
    [~, S_agg_next] = calculateVirtualSOC_agg(EVs, dt_long_minutes);

    %% 内层循环（短时间步长）
    for short_idx = 1:num_short_per_long
        step_idx = (long_idx - 1) * num_short_per_long + short_idx;
        t_relative_minute = t_long_start_minute + (short_idx - 1) * dt_minutes; 
    
        t_current_minute_abs = t_relative_minute + simulation_start_hour * 60; 
    
        current_absolute_hour = time_points_absolute(step_idx);

        temp_m3 = zeros(num_evs, 1);
        temp_S_original = zeros(num_evs, 1);
        temp_S_mod = zeros(num_evs, 1);
        temp_S_prime = zeros(num_evs, 1);
        temp_P_current = zeros(num_evs, 1);
        temp_E_current = zeros(num_evs, 1); % [新增] 临时存储当前实际电量
        temp_substate = zeros(num_evs, 1);
        
        % (新增) 为单体求和潜力创建临时存储
        temp_delta_p_plus_individual = zeros(num_evs, 1);
        temp_delta_p_minus_individual = zeros(num_evs, 1);

        %% 更新EV状态 (并行处理)
        EVs_in_parfor = EVs;

        parfor i = 1:num_evs
            EV = EVs_in_parfor(i);
            EV = updateLockState(EV, t_current_minute_abs);
            EV_for_handle = EV;
            EV_temp_with_handle = generateDemandCurve(EV_for_handle);
            current_P_val = 0;
            if isfield(EV_temp_with_handle, 'demandCurve') && isa(EV_temp_with_handle.demandCurve, 'function_handle')
                current_P_val = EV_temp_with_handle.demandCurve(lambda_star);
            else
                 switch EV.state
                    case 'LockON'
                        current_P_val = EV.P_N;
                    case {'LockOFF', 'OFF', 'ON'}
                        current_P_val = 0;
                    otherwise
                        current_P_val = 0;
                 end
            end
            EV.P_current = current_P_val;
            EV = updateLockState(EV, t_current_minute_abs);
            EV = calculateVirtualSOC_upgrade(EV, t_current_minute_abs, dt_minutes);

            if EV.t_dep > EV.t_in
                 m3_val = (EV.E_tar - EV.E_ini) / (EV.eta * ((EV.t_dep - EV.t_in) / 60));
            else
                 m3_val = 0;
            end
            temp_m3(i) = m3_val;

            temp_S_original(i) = EV.S_original;
            temp_S_mod(i) = EV.S_modified;
            temp_P_current(i) = EV.P_current;
            temp_E_current(i) = EV.E_actual; % [新增] 记录实际电量
            temp_S_prime(i) = getSPrime(EV);
            temp_substate(i) = strcmp(EV.substate, 'ON');
            % (新增) 计算单体调节潜力
            is_online_h = (current_absolute_hour >= (EV.t_in / 60)) && (current_absolute_hour < (EV.t_dep / 60)); % 判断是否在线 (小时)
            
            if EV.ptcp && is_online_h % 只有参与且在线的 EV 才有潜力
                 P_base_i = EV.P_base_sequence(step_idx); % 获取当前短步的基线功率
                 t_dep_h = EV.t_dep / 60; % 离网时间 (小时)

                 [DeltaP_plus_i, DeltaP_minus_i] = calculateEVAdjustmentPotentia_new(...
                     EV.E_reg_min, EV.E_reg_max, EV.E_actual, ... % 使用 E_actual 作为当前电量近似
                     t_dep_h, current_absolute_hour, ...
                     EV.P_N, P_base_i, EV.eta, t_adj); % 使用小时单位
                 
                 temp_delta_p_plus_individual(i) = DeltaP_plus_i;
                 temp_delta_p_minus_individual(i) = DeltaP_minus_i;
            else
                 % 未参与或不在线的EV，单体潜力为0
                 temp_delta_p_plus_individual(i) = 0;
                 temp_delta_p_minus_individual(i) = 0;
            end

            EVs_in_parfor(i) = EV;
        end % 结束 parfor i

        EVs = EVs_in_parfor; % 更新主 EVs 数组
        
        % (新增) 累加单体潜力
        individual_sum_DeltaP_plus = sum(temp_delta_p_plus_individual);
        individual_sum_DeltaP_minus = sum(temp_delta_p_minus_individual);

        %% (修改) 计算聚合潜力 (聚合模型) - 改进版：分组聚合
        active_participating_indices = find(arrayfun(@(ev, t_abs_h) ...
            ev.ptcp && (t_abs_h >= (ev.t_in / 60)) && (t_abs_h < (ev.t_dep / 60)), ...
            EVs, repmat(current_absolute_hour, num_evs, 1)));

        agg_DeltaP_plus = 0;
        agg_DeltaP_minus = 0;
        agg_P_real = 0; % [新增] 聚合模型实时功率

        if ~isempty(active_participating_indices)
            all_t_dep_h = [EVs(active_participating_indices).t_dep] / 60;
            t_rem_all = all_t_dep_h - current_absolute_hour;
            
            group_edges = [0, 1, 2, 4, 100]; 
            
            for g = 1:length(group_edges)-1
                group_mask = (t_rem_all >= group_edges(g)) & (t_rem_all < group_edges(g+1));
                group_indices_in_agg = active_participating_indices(group_mask);
                
                if ~isempty(group_indices_in_agg)
                    group_EVs = EVs(group_indices_in_agg);
                    
                    E_reg_min_agg = sum([group_EVs.E_reg_min]);
                    E_reg_max_agg = sum([group_EVs.E_reg_max]);
                    E_current_agg = sum([group_EVs.E_actual]); 
                    
                    t_dep_agg_h = mean([group_EVs.t_dep]) / 60; 
                    
                    p_on_agg = sum([group_EVs.P_N]);           
                    P_base_agg = sum(arrayfun(@(ev) ev.P_base_sequence(step_idx), group_EVs)); 
                    eta_agg = mean([group_EVs.eta]);           
        
                    [d_plus, d_minus] = calculateEVAdjustmentPotentia_new(...
                        E_reg_min_agg, E_reg_max_agg, E_current_agg, ...
                        t_dep_agg_h, current_absolute_hour, ...
                        p_on_agg, P_base_agg, eta_agg, t_adj);
                    
                    agg_DeltaP_plus = agg_DeltaP_plus + d_plus;
                    agg_DeltaP_minus = agg_DeltaP_minus + d_minus;
                    
                    % [新增] 计算该分组的实时功率求和，并累加到聚合模型实时功率
                    p_real_agg = sum([group_EVs.P_current]);
                    agg_P_real = agg_P_real + p_real_agg;
                end
            end
        end

        %% [新增] 计算特定组的聚合 SOC
        % 1. 提取该组 (PN=6.6, C=26.9, ptcp=1)
        current_group_EVs = EVs(group_indices);
        
        % 2. 筛选出当前处于“在线”状态的车辆 (LockON / ON)
        group_active_mask = strcmp({current_group_EVs.state}, 'ON') | strcmp({current_group_EVs.state}, 'LockON');
        active_group_EVs = current_group_EVs(group_active_mask);
        
        if ~isempty(active_group_EVs)
            sum_E_diff = 0;
            sum_Cr = 0;
            sum_S_prime = 0;
            
            for k = 1:length(active_group_EVs)
                ev_k = active_group_EVs(k);
                % 累加偏差能量 (备用)
                sum_E_diff = sum_E_diff + (ev_k.E_actual - ev_k.E_exp);
                sum_Cr = sum_Cr + (ev_k.C * mean(ev_k.r)); 
                
                % [关键修改] 累加 S' (直接使用 getSPrime 计算每辆车当前的 S')
                % 这样可以保证与 Lambda (S'域) 的比较是公平的，且消除了从 S_raw 反推 S' 的误差
                sum_S_prime = sum_S_prime + getSPrime(ev_k);
            end
            
            % 计算该组的加权平均原始 SOC
            results.S_agg_group_raw(step_idx) = - sum_E_diff / sum_Cr;
            
            % [关键修改] 直接使用平均 S' 作为输出
            results.S_agg_group_prime(step_idx) = sum_S_prime / length(active_group_EVs);
            
        else
            % 如果该组当前没有车在线
            results.S_agg_group_raw(step_idx) = NaN;
            results.S_agg_group_prime(step_idx) = NaN;
        end

        %% 记录结果
        results.EV_S_original(:, step_idx) = temp_S_original;
        results.EV_S_mod(:, step_idx) = temp_S_mod;
        results.m3 = temp_m3; 

        results.lambda(step_idx) = lambda_star;
        results.S_agg(step_idx) = S_agg_current;
        results.P_agg(step_idx) = sum(temp_P_current);
        results.P_agg_ptcp(step_idx) = sum(temp_P_current(ptcp_vec)); % [新增] 计算参与聚合的EV总功率
        results.P_cu(step_idx) = temp_P_current(10);

        results.EV_Up(step_idx)   = agg_DeltaP_plus;
        results.EV_Down(step_idx) = agg_DeltaP_minus;
        results.EV_Power(step_idx)= agg_P_real; % [新增] 记录聚合模型实时功率
        
        results.EV_Up_Individual_Sum(step_idx)   = individual_sum_DeltaP_plus;
        results.EV_Down_Individual_Sum(step_idx) = individual_sum_DeltaP_minus;
        
        results.EV_Up_Individual(:, step_idx) = temp_delta_p_plus_individual;
        results.EV_Down_Individual(:, step_idx) = temp_delta_p_minus_individual;
        results.SOC_EV(:, step_idx) = temp_S_original;
        
        results.EV_E_actual(:, step_idx) = temp_E_current; % [新增] 保存当前步的实际电量
        results.EV_S_prime(:, step_idx) = temp_S_prime; 
        results.EV_substate(:, step_idx) = temp_substate; 

    end % 结束 short_idx

    %% 更新聚合SOC
    S_agg_current = S_agg_next;

end % 结束 long_idx

results.P_tar = repelem(P_tar, num_short_per_long);

%% 结果保存与可视化
outputFileName = 'main_potential_5min_1000_8am.mat'; 
fprintf('\n正在保存结果到 %s ...\n', outputFileName);
try
    save(outputFileName, 'results', '-v7.3');
    fprintf('结果保存成功。\n');
catch ME_save
    fprintf('*** 保存结果文件时出错: %s ***\n', ME_save.message);
end

%% [独立绘图模块] 组内微观-宏观响应全景图 (单体S' vs 聚合S' vs Lambda)
% 目的: 验证参数完全相同的单体在 Lambda 控制下的分化与跟随行为

% 1. 检查是否存在特定组索引 (group_indices)
if ~exist('group_indices', 'var') || isempty(group_indices)
    warning('未找到 group_indices 变量，请确保主程序已运行且筛选了特定组。正在尝试重新筛选...');
    target_PN = 6.6; target_C = 26.86;
    group_indices = find(abs([EVs.P_N] - target_PN) < 1e-3 & ...
                         abs([EVs.C] - target_C) < 1e-3 & ...
                         [EVs.ptcp] == 1);
end

if ~isempty(group_indices)
    fprintf('正在绘制 %d 辆同质化单体 S'' 轨迹与 Lambda 的对比图...\n', length(group_indices));
    
    % 创建画布
    figure('Name', 'Group Micro-Macro Response', 'Color', 'w', 'Position', [100, 100, 1000, 600]);
    hold on;
    
    % --- Layer 1: 绘制组内每一辆单体的 S' 轨迹 (灰色细线) ---
    % 注意：results.EV_S_mod 存储的就是 S' (S_modified)
    % 为了避免图例爆炸，只给第一条线加标签，其他的 HandleVisibility off
    
    % 提取该组所有单体的 S' 数据矩阵
    group_s_prime_data = results.EV_S_mod(group_indices, :);
    
    % 绘制所有单体轨迹 (淡灰色，透明度低一点)
    h_ind = plot(time_points_absolute, group_s_prime_data, ...
        'Color', [0.7, 0.7, 0.7, 0.5], ... % 灰色，半透明
        'LineWidth', 0.5); 
    
    % 调整图例：将所有单体线的 HandleVisibility 设为 off
    for k = 1:length(h_ind)
        h_ind(k).HandleVisibility = 'off';
    end
    % 随便画一条假的灰色线用于生成图例
    h_legend_ind = plot(nan, nan, 'Color', [0.6, 0.6, 0.6], 'LineWidth', 1, 'DisplayName', '单体期望SOC');
    
    
    % --- Layer 2: 绘制该组的聚合平均 S' (蓝色粗实线) ---
    h_agg = plot(time_points_absolute, results.S_agg_group_prime, ...
        'Color', 'b', ...
        'LineWidth', 2.5, ...
        'DisplayName', '聚合均值期望SOC');
        
    % --- Layer 3: 绘制系统指令 Lambda (红色粗虚线) ---
    h_lambda = plot(time_points_absolute, results.lambda, ...
        'Color', 'r', ...
        'LineStyle', '--', ...
        'LineWidth', 2.5, ...
        'DisplayName', 'lambda');
    
    hold off;
    
    % --- 美化与标注 ---
    % 设置关注的时间段 (如 21:00 - 次日 04:00)，也可以全时段
    % xlim([21, 28]); 
    xlim([simulation_start_hour, simulation_end_hour]);
    
    xlabel('时间 (小时)', 'FontSize', 12, 'FontName', 'SimHei');
    ylabel('归一化期望SOC', 'FontSize', 12, 'FontName', 'SimHei');
    grid on;
    box on;
    set(gca, 'FontSize', 12);
    
    % 图例
    legend([h_legend_ind, h_agg, h_lambda], 'Location', 'best', 'FontSize', 11, 'FontName', 'SimHei');
    
    % --- 保存图片 ---
    output_img_name = 'lamda_soc_同组对比.png';
    print(gcf, output_img_name, '-dpng', '-r600');
    
    % 如果是 Windows，额外保存 EMF 矢量图
    if ispc
        print(gcf, 'lamda_soc_同组对比.emf', '-dmeta');
    end
    
    fprintf('绘图完成，已保存为 %s\n', output_img_name);
    
else
    warning('未找到符合条件的组，无法绘图。');
end