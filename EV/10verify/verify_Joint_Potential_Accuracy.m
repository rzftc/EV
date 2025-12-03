%% verify_Joint_Potential_Accuracy.m
% 联合验证程序：整合 EV(DER1), AC(DER2), PV(DER3), P2G(DER4), IDC(DER5), Industry(DER6)
clc; clear; close all;

set(0, 'DefaultAxesFontName', 'Microsoft YaHei');
set(0, 'DefaultTextFontName', 'Microsoft YaHei');

%% 1. 文件路径定义
der_files = {
    'DER1.mat', 'EV';
    'DER2.mat', 'AC';
    'DER3.mat', 'PV';
    'DER4.mat', 'P2G';
    'DER5.mat', 'IDC';
    'DER6.mat', 'Industry'
};

fprintf('==================================================\n');
fprintf('      多资源协同潜力预测联合验证      \n');
fprintf('==================================================\n');

%% 2. 数据加载与预处理容器
% 初始化存储结构
data_store = struct();
min_len = inf;

% 循环加载所有资源文件
for i = 1:size(der_files, 1)
    filename = der_files{i, 1};
    label = der_files{i, 2};
    
    fprintf('正在加载 %s (%s) ... ', filename, label);
    
    if ~exist(filename, 'file')
        error('未找到数据文件: %s', filename);
    end
    
    tmp = load(filename);
    if ~isfield(tmp, 'results')
        error('%s 中未找到 "results" 结构体。', filename);
    end
    res = tmp.results;
    
    % 提取数据并统一为列向量
    if strcmp(label, 'EV') % DER1 (特殊字段名)
        up_model = res.EV_Up(:);
        down_model = res.EV_Down(:);
        up_true = res.EV_Up_Individual_Sum(:);
        down_true = res.EV_Down_Individual_Sum(:);
        
    elseif strcmp(label, 'AC') % DER2 (特殊字段名)
        if isfield(res, 'Agg_Model_Potential_Up_History')
            up_model = res.Agg_Model_Potential_Up_History(:);
            down_model = res.Agg_Model_Potential_Down_History(:);
        else
            error('DER2 缺少聚合模型数据');
        end
        if isfield(res, 'Agg_P_Potential_Up_History')
            up_true = res.Agg_P_Potential_Up_History(:);
            down_true = res.Agg_P_Potential_Down_History(:);
        else
            error('DER2 缺少单体累加数据');
        end
        
    else % DER3 - DER6 (根据您要求的逻辑：聚合值 = 真值)
        switch label
            case 'PV'
                up_val = res.PV_Up(:);
                down_val = res.PV_Down(:);
            case 'P2G'
                up_val = res.P2G_Up(:);
                down_val = res.P2G_Down(:);
            case 'IDC'
                up_val = res.IDC_Up(:);
                down_val = res.IDC_Down(:);
            case 'Industry'
                up_val = res.Ind_Up(:);
                down_val = res.Ind_Down(:);
        end
        
        % 赋值：真值和聚合值设为相同
        up_model = up_val;
        up_true = up_val;
        down_model = down_val;
        down_true = down_val;
    end
    
    % 更新最小长度用于对齐
    min_len = min(min_len, length(up_model));
    
    % 存储到结构体
    data_store.(label).Up_Model = up_model;
    data_store.(label).Down_Model = down_model;
    data_store.(label).Up_True = up_true;
    data_store.(label).Down_True = down_true;
    
    fprintf('成功 (长度: %d)\n', length(up_model));
end

%% 3. 数据对齐与总量计算
fprintf('\n正在对齐数据 (截取前 %d 个时间步)...\n', min_len);

% 生成时间轴 (假设所有文件步长一致，使用5分钟步长)
dt_short = 5; 
time_axis = (0:min_len-1)' * dt_short / 60;

% 初始化总量
Total_Up_Model = zeros(min_len, 1);
Total_Down_Model = zeros(min_len, 1);
Total_Up_True = zeros(min_len, 1);
Total_Down_True = zeros(min_len, 1);

resources = fieldnames(data_store);
for i = 1:length(resources)
    label = resources{i};
    
    % 截取数据
    u_m = data_store.(label).Up_Model(1:min_len);
    d_m = data_store.(label).Down_Model(1:min_len);
    u_t = data_store.(label).Up_True(1:min_len);
    d_t = data_store.(label).Down_True(1:min_len);
    
    % 更新回结构体（用于绘图）
    data_store.(label).Up_Model = u_m;
    data_store.(label).Down_Model = d_m;
    data_store.(label).Up_True = u_t;
    data_store.(label).Down_True = d_t;
    
    % 累加总量
    Total_Up_Model = Total_Up_Model + u_m;
    Total_Down_Model = Total_Down_Model + d_m;
    Total_Up_True = Total_Up_True + u_t;
    Total_Down_True = Total_Down_True + d_t;
end

%% 4. 计算误差指标
% 计算比值 (真值 / 聚合值)
epsilon = 1e-3;

Ratio_Up = ones(size(Total_Up_Model));
valid_idx_up = abs(Total_Up_Model) > epsilon;
Ratio_Up(valid_idx_up) = Total_Up_True(valid_idx_up) ./ Total_Up_Model(valid_idx_up);

Ratio_Down = ones(size(Total_Down_Model));
valid_idx_down = abs(Total_Down_Model) > epsilon;
Ratio_Down(valid_idx_down) = Total_Down_True(valid_idx_down) ./ Total_Down_Model(valid_idx_down);

%% 5. 绘图展示
fprintf('正在生成对比图...\n');

% --- 图1: 六类资源的上调潜力堆叠图 (模型值) ---
figure('Color', 'w', 'Position', [100, 100, 1000, 600]);
hold on;
area_data_up = [];
legend_str = {};
for i = 1:length(resources)
    area_data_up = [area_data_up, data_store.(resources{i}).Up_Model];
    legend_str{end+1} = resources{i};
end
area(time_axis, area_data_up);
title('各类资源上调潜力构成 (聚合模型值)');
xlabel('时间 (小时)'); ylabel('功率 (kW)');
legend(legend_str, 'Location', 'best');
grid on; set(gca, 'FontSize', 12);

% --- 图2: 六类资源的下调潜力堆叠图 (模型值) ---
figure('Color', 'w', 'Position', [150, 150, 1000, 600]);
hold on;
area_data_down = [];
for i = 1:length(resources)
    area_data_down = [area_data_down, data_store.(resources{i}).Down_Model];
end
area(time_axis, area_data_down);
title('各类资源下调潜力构成 (聚合模型值)');
xlabel('时间 (小时)'); ylabel('功率 (kW)');
legend(legend_str, 'Location', 'best');
grid on; set(gca, 'FontSize', 12);

% --- 图3: 联合系统总量对比 (上调) ---
figure('Color', 'w', 'Position', [200, 200, 800, 500]);
plot(time_axis, Total_Up_Model, 'r-', 'LineWidth', 1.5, 'DisplayName', '总聚合模型预测');
hold on;
plot(time_axis, Total_Up_True, 'k--', 'LineWidth', 1.5, 'DisplayName', '总单体累加真值');
title('联合系统上调潜力对比');
xlabel('时间 (小时)'); ylabel('功率 (kW)');
legend('Location', 'best'); grid on; set(gca, 'FontSize', 12);

% --- 图4: 联合系统总量对比 (下调) ---
figure('Color', 'w', 'Position', [250, 250, 800, 500]);
plot(time_axis, Total_Down_Model, 'b-', 'LineWidth', 1.5, 'DisplayName', '总聚合模型预测');
hold on;
plot(time_axis, Total_Down_True, 'k--', 'LineWidth', 1.5, 'DisplayName', '总单体累加真值');
title('联合系统下调潜力对比');
xlabel('时间 (小时)'); ylabel('功率 (kW)');
legend('Location', 'best'); grid on; set(gca, 'FontSize', 12);

% --- 图5: 准确度比值分析 ---
figure('Color', 'w', 'Position', [300, 300, 800, 500]);
subplot(2,1,1);
plot(time_axis, Ratio_Up, 'm-', 'LineWidth', 1.5);
yline(1.0, 'k--', 'LineWidth', 1.5);
title('上调潜力准确度比值 (真值/模型值)');
xlabel('时间 (小时)'); ylabel('比值'); ylim([0.5, 1.5]); grid on;

subplot(2,1,2);
plot(time_axis, Ratio_Down, 'c-', 'LineWidth', 1.5);
yline(1.0, 'k--', 'LineWidth', 1.5);
title('下调潜力准确度比值 (真值/模型值)');
xlabel('时间 (小时)'); ylabel('比值'); ylim([0.5, 1.5]); grid on;

%% 6. 总量偏差计算报告
fprintf('\n========== 分布式可控资源联合潜力预测验证报告 ==========\n');

% 上调总量
Sum_Model_Up = sum(Total_Up_Model);
Sum_True_Up  = sum(Total_Up_True);

% 下调总量
Sum_Model_Down = sum(Total_Down_Model);
Sum_True_Down  = sum(Total_Down_True);

% 计算相对误差
if abs(Sum_True_Up) > 1e-3
    Err_Up = abs(Sum_Model_Up - Sum_True_Up) / abs(Sum_True_Up) * 100;
else
    Err_Up = 0;
end

if abs(Sum_True_Down) > 1e-3
    Err_Down = abs(Sum_Model_Down - Sum_True_Down) / abs(Sum_True_Down) * 100;
else
    Err_Down = 0;
end

fprintf('【联合上调潜力】\n');
fprintf('  - 聚合模型总量:                                %.2f kW·step\n', Sum_Model_Up);
fprintf('  - 单体累加总量:                                %.2f kW·step\n', Sum_True_Up);
fprintf('  - 可控资源与电网时空互济规划态潜力预测平均精确度  %.2f%%\n', 100-Err_Up);

fprintf('【联合下调潜力】\n');
fprintf('  - 聚合模型总量:                                %.2f kW·step\n', Sum_Model_Down);
fprintf('  - 单体累加总量:                                %.2f kW·step\n', Sum_True_Down);
fprintf('  - 可控资源与电网时空互济规划态潜力预测平均精确度: %.2f%%\n', 100-Err_Down);

fprintf('======================================================\n');