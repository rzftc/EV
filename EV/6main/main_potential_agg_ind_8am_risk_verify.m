%% main_potential_agg_ind_8am_risk_verify.m
% (修改版: 增加EV用户失约/违约行为的显性验证)
% 1. 生成用户失约场景 (提前离网)
% 2. 区分“聚合商视角(合约)”与“物理视角(实际)”的潜力计算
% 3. 可视化违约风险造成的容量缺口
% 4. [新增] 可视化违约风险造成的聚合功率偏差，并保存高DPI图片

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
    'P_agg_actual',  zeros(1, total_steps), ... % [新增] 考虑失约后的实际聚合功率
    'P_base',        zeros(1, total_steps), ...
    'S_agg',         zeros(1, total_steps), ...
    'lambda',        zeros(1, total_steps), ...
    'EV_S_original', zeros(num_evs, total_steps), ...
    'EV_S_mod',      zeros(num_evs, total_steps), ...
    'm3',            zeros(num_evs, 1), ...
    'P_tar',         zeros(1, total_steps), ...
    'P_cu',          zeros(1, total_steps), ...
    'EV_Up',         zeros(1, total_steps), ... % 聚合模型上调潜力 (聚合商视角)
    'EV_Down',       zeros(1, total_steps), ... % 聚合模型下调潜力 (聚合商视角)
    'EV_Power',      zeros(1, total_steps), ... % 聚合模型实时功率
    'EV_Up_Individual_Sum',   zeros(1, total_steps), ... % 单体求和上调潜力 (物理视角)
    'EV_Down_Individual_Sum', zeros(1, total_steps),  ... % 单体求和下调潜力 (物理视角)
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

%% [关键新增] 生成用户失约/违约场景数据 (基于论文1.3.3节)
fprintf('正在生成用户失约/违约场景数据...\n');
% 定义失约风险参数
risk_early_departure_prob = 0.3; % 30%的用户可能会提前离网(失约)
risk_mu = 60;    % 平均提前离网时间 (分钟)
risk_sigma = 20; % 标准差 (分钟)

% 初始化合约时间（即聚合商以为的时间，等于原计划时间）
temp_t_dep_contract = [EVs.t_dep]'; 
temp_t_in = [EVs.t_in]';

% 初始化实际时间（物理上真实发生的时间）
temp_t_dep_actual = temp_t_dep_contract;

% 生成随机失约行为
is_breach_user = rand(num_evs, 1) < risk_early_departure_prob; % 标记违约用户
early_minutes = max(0, normrnd(risk_mu, risk_sigma, num_evs, 1)); % 生成提前时间

% 更新违约用户的实际离网时间 (实际时间 = 合约时间 - 提前时间)
temp_t_dep_actual(is_breach_user) = temp_t_dep_contract(is_breach_user) - early_minutes(is_breach_user);

% 确保实际离网时间不早于入网时间 (至少停15分钟)
temp_t_dep_actual = max(temp_t_dep_actual, temp_t_in + 15);

% 将数据存回EVs结构体，方便后续调用
for i = 1:num_evs
    EVs(i).t_dep_contract = temp_t_dep_contract(i); % 合约时间 (用于聚合计算)
    EVs(i).t_dep_actual   = temp_t_dep_actual(i);   % 实际时间 (用于单体物理更新)
end
fprintf('已生成失约场景：%.1f%% 的用户提前离网，平均提前 %.1f 分钟。\n', ...
    mean(is_breach_user)*100, mean(early_minutes(is_breach_user)));


%% 外层循环（长时间步长）
for long_idx = 1:num_long_steps
    t_long_start_minute = (long_idx - 1) * dt_long_minutes;

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
        temp_P_current = zeros(num_evs, 1);
        temp_P_actual = zeros(num_evs, 1); % [新增] 临时存储当前实际物理功率(考虑失约)
        temp_E_current = zeros(num_evs, 1); % [新增] 临时存储当前实际电量
        
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

            EV = calculateVirtualSOC_upgrade(EV, t_current_minute_abs, dt_minutes);

            if EV.t_dep > EV.t_in
                 m3_val = (EV.E_tar - EV.E_ini) / (EV.eta * ((EV.t_dep - EV.t_in) / 60));
            else
                 m3_val = 0;
            end
            temp_m3(i) = m3_val;

            temp_S_original(i) = EV.S_original;
            temp_S_mod(i) = EV.S_modified;
            temp_P_current(i) = EV.P_current; % 合约视角下的功率
            temp_E_current(i) = EV.E_actual; % [新增] 记录实际电量
            
            % [修改] 计算单体调节潜力 (物理视角 - 使用实际离网时间)
            t_dep_actual_h = EV.t_dep_actual / 60;
            is_physically_online = (current_absolute_hour >= (EV.t_in / 60)) && ...
                                   (current_absolute_hour < t_dep_actual_h);
            
            % [新增] 计算实际物理功率 (如果已经提前离网，则功率为0)
            if is_physically_online
                temp_P_actual(i) = EV.P_current;
            else
                temp_P_actual(i) = 0;
            end
            
            if EV.ptcp && is_physically_online 
                 P_base_i = EV.P_base_sequence(step_idx); % 获取当前短步的基线功率
                 
                 % [关键修改] 使用 t_dep_actual_h 进行潜力计算，反映物理真实情况
                 [DeltaP_plus_i, DeltaP_minus_i] = calculateEVAdjustmentPotentia_new(...
                     EV.E_reg_min, EV.E_reg_max, EV.E_actual, ... 
                     t_dep_actual_h, current_absolute_hour, ...
                     EV.P_N, P_base_i, EV.eta, t_adj); 
                 
                 temp_delta_p_plus_individual(i) = DeltaP_plus_i;
                 temp_delta_p_minus_individual(i) = DeltaP_minus_i;
            else
                 % 物理上已离网或不参与，潜力为0
                 temp_delta_p_plus_individual(i) = 0;
                 temp_delta_p_minus_individual(i) = 0;
            end

            EVs_in_parfor(i) = EV;
        end % 结束 parfor i

        EVs = EVs_in_parfor; % 更新主 EVs 数组
        
        % (新增) 累加单体潜力 (这是物理真实的潜力)
        individual_sum_DeltaP_plus = sum(temp_delta_p_plus_individual);
        individual_sum_DeltaP_minus = sum(temp_delta_p_minus_individual);

        %% (修改) 计算聚合潜力 (聚合模型) - 改进版：分组聚合
        % [关键修改] 聚合商使用合约时间 t_dep_contract 进行筛选
        active_participating_indices = find(arrayfun(@(ev, t_abs_h) ...
            ev.ptcp && (t_abs_h >= (ev.t_in / 60)) && (t_abs_h < (ev.t_dep_contract / 60)), ...
            EVs, repmat(current_absolute_hour, num_evs, 1)));

        agg_DeltaP_plus = 0;
        agg_DeltaP_minus = 0;
        agg_P_real = 0; % [新增] 聚合模型实时功率

        if ~isempty(active_participating_indices)
            % [关键修改] 聚合计算时使用 t_dep_contract (聚合商以为的时间)
            all_t_dep_h = [EVs(active_participating_indices).t_dep_contract] / 60;
            t_rem_all = all_t_dep_h - current_absolute_hour;
            
            group_edges = [0, 1, 2, 4, 100]; 
            
            for g = 1:length(group_edges)-1
                group_mask = (t_rem_all >= group_edges(g)) & (t_rem_all < group_edges(g+1));
                group_indices = active_participating_indices(group_mask);
                
                if ~isempty(group_indices)
                    group_EVs = EVs(group_indices);
                    
                    E_reg_min_agg = sum([group_EVs.E_reg_min]);
                    E_reg_max_agg = sum([group_EVs.E_reg_max]);
                    E_current_agg = sum([group_EVs.E_actual]); 
                    
                    % [关键修改] 使用合约时间计算聚合离网时间
                    t_dep_agg_h = mean([group_EVs.t_dep_contract]) / 60; 
                    
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

        %% 记录结果
        results.EV_S_original(:, step_idx) = temp_S_original;
        results.EV_S_mod(:, step_idx) = temp_S_mod;
        results.m3 = temp_m3; 

        results.lambda(step_idx) = lambda_star;
        results.S_agg(step_idx) = S_agg_current;
        results.P_agg(step_idx) = sum(temp_P_current); % 合约总功率
        results.P_agg_actual(step_idx) = sum(temp_P_actual); % [新增] 实际总功率(考虑失约)
        results.P_agg_ptcp(step_idx) = sum(temp_P_current(ptcp_vec)); % [新增] 计算参与聚合的EV总功率
        results.P_cu(step_idx) = temp_P_current(10);

        results.EV_Up(step_idx)   = agg_DeltaP_plus;  % 聚合模型估算值 (高估)
        results.EV_Down(step_idx) = agg_DeltaP_minus; % 聚合模型估算值 (高估)
        results.EV_Power(step_idx)= agg_P_real; 
        
        results.EV_Up_Individual_Sum(step_idx)   = individual_sum_DeltaP_plus;  % 物理真实值 (计及违约)
        results.EV_Down_Individual_Sum(step_idx) = individual_sum_DeltaP_minus; % 物理真实值 (计及违约)
        
        results.EV_Up_Individual(:, step_idx) = temp_delta_p_plus_individual;
        results.EV_Down_Individual(:, step_idx) = temp_delta_p_minus_individual;
        results.SOC_EV(:, step_idx) = temp_S_original;
        
        results.EV_E_actual(:, step_idx) = temp_E_current; % [新增] 保存当前步的实际电量

    end % 结束 short_idx

    %% 更新聚合SOC
    S_agg_current = S_agg_next;

end % 结束 long_idx

results.P_tar = repelem(P_tar, num_short_per_long);

%% 结果保存与可视化
outputFileName = 'main_potential_5min_1000_8am_risk_verify.mat'; 
fprintf('\n正在保存结果到 %s ...\n', outputFileName);
try
    save(outputFileName, 'results', '-v7.3');
    fprintf('结果保存成功。\n');
catch ME_save
    fprintf('*** 保存结果文件时出错: %s ***\n', ME_save.message);
end

% 设置绘图字体
defaultFont = 'Microsoft YaHei';

% 1. 可视化失约风险验证 (调节能力对比)
f1 = figure('Name', 'EV用户失约行为对调节能力的影响验证', 'Color', 'w');
% 上调潜力对比
subplot(2,1,1);
plot(time_points_absolute, results.EV_Up, 'r--', 'LineWidth', 1.5, 'DisplayName', '合约潜力');
hold on;
plot(time_points_absolute, results.EV_Up_Individual_Sum, 'r-', 'LineWidth', 1.5, 'DisplayName', '失约潜力');
% 填充差值区域 (风险)
x_conf = [time_points_absolute, fliplr(time_points_absolute)];
y_conf = [results.EV_Up, fliplr(results.EV_Up_Individual_Sum)];
fill(x_conf, y_conf, 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'DisplayName', '失约造成的容量缺口');
ylabel('上调潜力 (kW)', 'FontName', defaultFont);
title('用户失约行为下的调节能力偏差验证 (上调)', 'FontName', defaultFont);
legend('Location', 'best', 'FontName', defaultFont); grid on; 
xlim([simulation_start_hour, simulation_end_hour]);
set(gca, 'FontName', defaultFont);

% 下调潜力对比
subplot(2,1,2);
plot(time_points_absolute, results.EV_Down, 'b--', 'LineWidth', 1.5, 'DisplayName', '合约潜力');
hold on;
plot(time_points_absolute, results.EV_Down_Individual_Sum, 'b-', 'LineWidth', 1.5, 'DisplayName', '失约潜力');
% 填充差值区域
y_conf_down = [results.EV_Down, fliplr(results.EV_Down_Individual_Sum)];
fill(x_conf, y_conf_down, 'b', 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'DisplayName', '失约造成的容量缺口');
ylabel('下调潜力 (kW)', 'FontName', defaultFont); xlabel('时间 (小时)', 'FontName', defaultFont);
title('用户失约行为下的调节能力偏差验证 (下调)', 'FontName', defaultFont);
legend('Location', 'best', 'FontName', defaultFont); grid on; 
xlim([simulation_start_hour, simulation_end_hour]);
set(gca, 'FontName', defaultFont);

% 保存图1
exportgraphics(f1, 'EV用户失约调节能力验证.png', 'Resolution', 300);
exportgraphics(f1, 'EV用户失约调节能力验证.emf', 'Resolution', 300);

% 2. [新增] 可视化失约对聚合功率的影响
f2 = figure('Name', 'EV用户失约行为对聚合功率的影响验证', 'Color', 'w');

% 1) 绘制合约预估功率 (蓝色虚线 - 代表理想情况)
plot(time_points_absolute, results.P_agg, 'b--', 'LineWidth', 1.5, ...
    'DisplayName', '合约功率');
hold on;

% 2) 绘制实际聚合功率 (红色实线 - 代表风险情况)
plot(time_points_absolute, results.P_agg_actual, 'r-', 'LineWidth', 1.5, ...
    'DisplayName', '失约功率');

% 3) 填充差值区域 (红色半透明 - 突出缺口)
% 构造填充多边形的顶点
x_conf = [time_points_absolute, fliplr(time_points_absolute)];
y_conf_power = [results.P_agg, fliplr(results.P_agg_actual)];
fill(x_conf, y_conf_power, 'r', 'FaceAlpha', 0.15, 'EdgeColor', 'none', ...
    'DisplayName', '失约造成的功率偏差');

% 4) 坐标轴设置
ylabel('聚合功率 (kW)', 'FontName', defaultFont, 'FontSize', 12); 
xlabel('时间 (小时)', 'FontName', defaultFont, 'FontSize', 12);
title('用户失约行为下的聚合功率偏差验证', 'FontName', defaultFont, 'FontSize', 14);
legend('Location', 'best', 'FontName', defaultFont, 'FontSize', 10); 
grid on; box on;

% 关键修改：强制横坐标范围从 8:00 (Day 1) 到 32:00 (Day 2, 即次日8:00)
xlim([8, 32]); 
% 设置刻度：每4小时一个刻度，更清晰
xticks(8:4:32);
xticklabels({'08:00','12:00','16:00','20:00','00:00','04:00','08:00'});

set(gca, 'FontName', defaultFont, 'FontSize', 11);

% 保存图2
exportgraphics(f2, 'EV用户失约聚合功率验证.png', 'Resolution', 300);
exportgraphics(f2, 'EV用户失约聚合功率验证.emf', 'Resolution', 300);