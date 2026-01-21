%% main_potential_agg_ind_all.m
% (修改版: 批量仿真不同激励价格 0-50，并分别保存结果)
% (包含: 单体求和潜力(旧版算法) vs 聚合模型潜力，分组聚合改进)
% (本次修正: 恢复基于偏差度 r 的原版计算理论，参数随机生成)

clc; clear; close all;

%% ===================== 批量仿真设置 =====================
% 定义10个均匀分布的激励价格 (0 到 50)
incentive_prices = linspace(0, 50, 10);

% 开始批量仿真循环
for sim_idx = 1:length(incentive_prices)
    current_incentive_price = incentive_prices(sim_idx);
    rng(2024);

    fprintf('\n==========================================================\n');
    fprintf('正在执行第 %d/%d 次仿真: 统一激励价格 = %.2f\n', sim_idx, length(incentive_prices), current_incentive_price);
    fprintf('==========================================================\n');

    %% 初始化参数
    excelFile = 'evtest_8am.xlsx';
    if ~exist(excelFile, 'file')
        generateEVParameters_real(excelFile, 1000, 1.0);
        fprintf('已生成参数模板: %s\n', excelFile);
    end
    % 每次循环重新加载，确保EV结构体状态重置
    [EVs, t_sim, ~, ~, P_tar] = initializeFromExcel_8am(excelFile);
    fprintf('成功加载%d辆EV数据\n', length(EVs));

    %% 时间参数定义
    dt_short = 5;     % 短时间步长 (分钟)
    dt_long = 15;       % 长时间步长 (分钟)
    simulation_start_hour = 8;
    simulation_end_hour   = 32;
    dt = dt_short / 60;       % 短步长 (小时)
    dt_minutes = dt_short;    % 短步长 (分钟)
    dt_long_minutes = dt_long;% 长步长 (分钟)
    t_adj = 5 / 60;           % [修正] 调节时长改为 5分钟 (与bound版保持一致以匹配旧算法)

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
        'EV_Up_Individual_Sum',   zeros(1, total_steps), ... % 单体求和上调潜力
        'EV_Down_Individual_Sum', zeros(1, total_steps),  ... % 单体求和下调潜力
        ...
        'SOC_EV',             zeros(num_evs, total_steps), ...
        'EV_Up_Individual',   zeros(num_evs, total_steps), ...
        'EV_Down_Individual', zeros(num_evs, total_steps), ...
        ...
        'EV_E_actual',        zeros(num_evs, total_steps), ...
        'EV_E_baseline',      zeros(num_evs, total_steps), ...
        'P_base_agg',         zeros(1, total_steps), ...
        'EV_t_in',            zeros(num_evs, 1), ...
        'EV_t_dep',           zeros(num_evs, 1) ...
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
    r_vec = [EVs.r]'; % [新增] 获取初始 r

    E_tar_original_vec = E_tar_set_vec;
    SOC_original_vec = zeros(num_evs, 1);
    E_reg_min_vec = E_tar_original_vec;
    E_reg_max_vec = E_tar_original_vec;

    p_min = 15; p_max = 50;
    p_min_prime = 10; p_max_prime = 40;
    base_Price_vec = 26 * ones(num_evs, 1);

    % ---------------------------------------------------------
    % 使用当前循环的激励价格覆盖原有数据
    % ---------------------------------------------------------
    p_incentive_vec = current_incentive_price * ones(num_evs, 1);

    % 将统一后的价格写回结构体（保持数据一致性）
    p_inc_cell = num2cell(p_incentive_vec);
    [EVs.p_incentive] = p_inc_cell{:};
    % ---------------------------------------------------------

    % 基于新价格计算参与率和灵活性
    participation_probabilities_vec = calculateParticipation(p_incentive_vec, base_Price_vec);
    ptcp_vec = (rand(num_evs, 1) < participation_probabilities_vec);

    E_tar_max_flex_vec = 0.2 * C_vec;
    % 计算理论最大偏差范围
    [deltaE_up_vec, deltaE_down_vec] = incentiveTempEV_updown(p_incentive_vec, p_min, p_max, p_min_prime, p_max_prime, E_tar_max_flex_vec);

    participating_indices = find(ptcp_vec);
    if ~isempty(participating_indices)
        % =========================================================
        % [修改点 - 参数逻辑] deltaE 在 deltaE_down 和 deltaE_up 之间随机生成
        % 并修正 r 值
        % =========================================================
        min_safe_r = 0.05;
        idx = participating_indices;
        d_up = deltaE_up_vec(idx);
        d_down = deltaE_down_vec(idx);
        
        % 生成随机实际偏差: actual = down + (up - down) * rand
        actual_delta_E = d_down + (d_up - d_down) .* rand(length(idx), 1);
        
        % 更新 E_reg_min 和 E_reg_max (对称应用)
        E_reg_min_vec(idx) = E_tar_original_vec(idx) - actual_delta_E;
        E_reg_max_vec(idx) = E_tar_original_vec(idx) + actual_delta_E;
        
        % 物理边界截断
        E_reg_min_vec(idx) = max(E_reg_min_vec(idx), E_ini_vec(idx));
        E_reg_max_vec(idx) = min(E_reg_max_vec(idx), C_vec(idx));
        
        % 修正 r: r = actual_delta_E / C
        r_vec(idx) = max(actual_delta_E ./ C_vec(idx), min_safe_r);
        
        % 存回 deltaE 供参考
        deltaE_up_vec(idx) = actual_delta_E;
        deltaE_down_vec(idx) = actual_delta_E;
    end

    % 写回EVs结构体
    E_tar_orig_cell = num2cell(E_tar_original_vec);
    SOC_orig_cell = num2cell(SOC_original_vec);
    ptcp_cell = num2cell(ptcp_vec);
    E_reg_min_cell = num2cell(E_reg_min_vec);
    E_reg_max_cell = num2cell(E_reg_max_vec);
    r_cell = num2cell(r_vec); % [新增] 准备写回修正后的 r

    [EVs.E_tar_original] = E_tar_orig_cell{:};
    [EVs.SOC_original] = SOC_orig_cell{:};
    [EVs.ptcp] = ptcp_cell{:};
    [EVs.E_reg_min] = E_reg_min_cell{:};
    [EVs.E_reg_max] = E_reg_max_cell{:};
    [EVs.r] = r_cell{:}; % [新增] 更新结构体中的 r (这对旧版算法至关重要)

    % [新增] 将 deltaE 参数写入 EVs (供参考)
    deltaE_up_cell = num2cell(deltaE_up_vec);
    deltaE_down_cell = num2cell(deltaE_down_vec);
    [EVs.deltaE_up] = deltaE_up_cell{:};
    [EVs.deltaE_down] = deltaE_down_cell{:};

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
    [EVs.r] = r_cell{:}; % 再次确认 r 已更新

    %% 外层循环（长时间步长）
    for long_idx = 1:num_long_steps
        t_long_start_minute = (long_idx - 1) * dt_long_minutes;

        %% 长时间步处理
        [lambda_star] = aggregateEVs(EVs, P_tar(long_idx));
        [~, S_agg_next] = calculateVirtualSOC_agg(EVs, dt_long_minutes);

        %% 内层循环（短时间步长）
        for short_idx = 1:num_short_per_long
            step_idx = (long_idx - 1) * num_short_per_long + short_idx;
            
            % --- [修改点 1] 时间计算修正 ---
            t_relative_minute = t_long_start_minute + (short_idx - 1) * dt_minutes;
            t_current_minute_abs = t_relative_minute + simulation_start_hour * 60;
            current_absolute_hour = time_points_absolute(step_idx);
            % ---------------------------

            temp_m3 = zeros(num_evs, 1);
            temp_S_original = zeros(num_evs, 1);
            temp_S_mod = zeros(num_evs, 1);
            temp_P_current = zeros(num_evs, 1);
            temp_E_current = zeros(num_evs, 1);

            % 为单体求和潜力创建临时存储
            temp_delta_p_plus_individual = zeros(num_evs, 1);
            temp_delta_p_minus_individual = zeros(num_evs, 1);

            %% 更新EV状态 (并行处理)
            EVs_in_parfor = EVs;

            parfor i = 1:num_evs
                EV = EVs_in_parfor(i);
                if EV.ptcp
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
                    temp_P_current(i) = EV.P_current;
                    temp_E_current(i) = EV.E_actual;

                    % 计算单体调节潜力
                    is_online_h = (current_absolute_hour >= (EV.t_in / 60)) && (current_absolute_hour < (EV.t_dep / 60));

                    if EV.ptcp && is_online_h
                        P_base_i = EV.P_base_sequence(step_idx);
                        % =========================================================
                        % [修改点 - 潜力计算] 使用旧版 calculateEVAdjustmentPotentia
                        % =========================================================
                        t_dep_h = EV.t_dep / 60;
                        t_in_h = EV.t_in / 60;
                        P_max_i = EV.P_N;
                        P_min_i = 0; % 单向充电
                        
                        [DeltaP_plus_i, DeltaP_minus_i] = calculateEVAdjustmentPotentia(...
                             EV.C, EV.r, EV.eta, EV.E_tar, EV.E_ini, EV.E_actual, ...
                             t_dep_h, t_in_h, P_max_i, P_min_i, P_base_i, ...
                             EV.S_modified, t_adj);

                        temp_delta_p_plus_individual(i) = DeltaP_plus_i;
                        temp_delta_p_minus_individual(i) = DeltaP_minus_i;
                    else
                        temp_delta_p_plus_individual(i) = 0;
                        temp_delta_p_minus_individual(i) = 0;
                    end

                    EVs_in_parfor(i) = EV;
                end
            end % 结束 parfor i

            EVs = EVs_in_parfor;

            % 累加单体潜力
            individual_sum_DeltaP_plus = sum(temp_delta_p_plus_individual);
            individual_sum_DeltaP_minus = sum(temp_delta_p_minus_individual);

            %% 计算聚合潜力 (聚合模型) - 分组聚合
            % (此处保持原有的聚合模型逻辑不变，仅单体部分按照要求修改为旧算法，
            % 这样可以对比"旧版单体累加"与"新版分组聚合"的差异)
            active_participating_indices = find(arrayfun(@(ev, t_abs_h) ...
                ev.ptcp && (t_abs_h >= (ev.t_in / 60)) && (t_abs_h < (ev.t_dep / 60)), ...
                EVs, repmat(current_absolute_hour, num_evs, 1)));

            agg_DeltaP_plus = 0;
            agg_DeltaP_minus = 0;

            if ~isempty(active_participating_indices)
                all_t_dep_h = [EVs(active_participating_indices).t_dep] / 60;
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
                    end
                end
            end

            %% 记录结果
            results.EV_S_original(:, step_idx) = temp_S_original;
            results.EV_S_mod(:, step_idx) = temp_S_mod;
            results.m3 = temp_m3;

            results.lambda(step_idx) = lambda_star;
            results.S_agg(step_idx) = S_agg_current;
            results.P_agg(step_idx) = sum(temp_P_current);
            results.P_cu(step_idx) = temp_P_current(10);

            results.EV_Up(step_idx)   = agg_DeltaP_plus;
            results.EV_Down(step_idx) = agg_DeltaP_minus;

            results.EV_Up_Individual_Sum(step_idx)   = individual_sum_DeltaP_plus;
            results.EV_Down_Individual_Sum(step_idx) = individual_sum_DeltaP_minus;

            results.EV_Up_Individual(:, step_idx) = temp_delta_p_plus_individual;
            results.EV_Down_Individual(:, step_idx) = temp_delta_p_minus_individual;
            results.SOC_EV(:, step_idx) = temp_S_original;

            results.EV_E_actual(:, step_idx) = temp_E_current;

        end % 结束 short_idx

        %% 更新聚合SOC
        S_agg_current = S_agg_next;

    end % 结束 long_idx

    results.P_tar = repelem(P_tar, num_short_per_long);

    %% 结果保存
    % 根据当前价格生成唯一文件名
    outputFileName = sprintf('results_incentive_%.2f_1000_bound.mat', current_incentive_price);
    fprintf('\n正在保存结果到 %s ...\n', outputFileName);
    try
        save(outputFileName, 'results', '-v7.3');
        fprintf('结果保存成功。\n');
    catch ME_save
        fprintf('*** 保存结果文件时出错: %s ***\n', ME_save.message);
    end

end % 结束 sim_idx 循环

fprintf('\n==========================================================\n');
fprintf('所有 10 组仿真已全部完成。\n');
fprintf('==========================================================\n');