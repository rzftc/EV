%% 工业大负荷时空互济潜力预测模型 (宏观场景版)
% 基于历史典型日曲线与未来全网负荷预测数据
% 数据来源：附件1 CSV文件 & 用户提供的2024/2025/2030规划数据
% 修改时间：2025-12-08

%% 1. 清理工作区
clear; clc; close all;

%% 2. 基础参数设置 (根据用户提供数据)
% 单位统一转换为：千瓦 (kW)
Unit_WanKW_to_KW = 10^4;
Unit_YiKW_to_KW = 10^8;

% 全网最大负荷预测 (亿千瓦 -> kW)
Grid_Max_Load_2024 = 1.145 * Unit_YiKW_to_KW;
Grid_Max_Load_2025 = 1.302 * Unit_YiKW_to_KW;
Grid_Max_Load_2030 = 1.600 * Unit_YiKW_to_KW;

% 关键比例参数
Industrial_Ratio_Min = 0.60;
Industrial_Ratio_Max = 0.70;
Industrial_Ratio_Avg = 0.65; % 仿真采用平均值 65%

Regulation_Rate = 0.05;      % 最大调节能力占比 5%

%% 3. 读取历史典型数据 (提取波形)
csv_filename = '附件1 工业负荷曲线.xlsx';

% 检查文件是否存在
if ~isfile(csv_filename)
    % 如果脚本在LargeLoad目录下，尝试向上寻找或直接指定绝对路径
    if isfile(['../' csv_filename])
        csv_filename = ['../' csv_filename];
    else
        error('未找到数据文件: %s，请确认文件路径。', csv_filename);
    end
end

% 读取数据
opts = detectImportOptions(csv_filename);
opts.VariableNamingRule = 'preserve';
data_table = readtable(csv_filename, opts);

% 提取24小时负荷数据 (第4列到第27列对应00:00-23:00)
% 注意：CSV中单位为"万千瓦"
raw_load_data_WanKW = table2array(data_table(:, 4:27));
raw_load_data_KW = raw_load_data_WanKW * Unit_WanKW_to_KW;

% 提取典型日波形：选取历史数据中负荷峰值最大的一天作为“典型极限场景”
[max_vals, day_indices] = max(max(raw_load_data_KW, [], 2));
typical_day_load_KW = raw_load_data_KW(day_indices, :);

% 归一化处理：获取单纯的波形系数 (0~1)
typical_profile_norm = typical_day_load_KW / max(typical_day_load_KW);

% 时间轴
t_axis = 0:23;

%% 4. 多场景仿真计算
scenarios = [
    2024, Grid_Max_Load_2024;
    2025, Grid_Max_Load_2025;
    2030, Grid_Max_Load_2030
];

results = struct();

fprintf('---------------- 仿真计算结果 ----------------\n');
for i = 1:size(scenarios, 1)
    year = scenarios(i, 1);
    grid_max = scenarios(i, 2);
    
    % 1. 计算该年份预测的工业负荷峰值
    ind_peak_load = grid_max * Industrial_Ratio_Avg;
    
    % 2. 生成该年份的工业负荷曲线 (波形 * 峰值)
    ind_load_curve = typical_profile_norm * ind_peak_load;
    
    % 3. 计算调节能力曲线 (负荷 * 5%)
    reg_capacity_curve = ind_load_curve * Regulation_Rate;
    
    % 存储结果
    results.(['y' num2str(year)]).load = ind_load_curve;
    results.(['y' num2str(year)]).reg = reg_capacity_curve;
    results.(['y' num2str(year)]).peak = ind_peak_load;
    results.(['y' num2str(year)]).max_reg = max(reg_capacity_curve);
    
    % 打印输出
    fprintf('[%d年] 全网最大负荷: %.3f 亿kW | 工业峰值: %.3f 亿kW | 最大调节潜力: %.3f 万kW\n', ...
        year, grid_max/Unit_YiKW_to_KW, ind_peak_load/Unit_YiKW_to_KW, max(reg_capacity_curve)/Unit_WanKW_to_KW);
end
fprintf('--------------------------------------------\n');

%% 5. 可视化绘图

% 设置绘图风格
figureConfig = {'Color','white', 'Position',[100 100 1400 600]};
fontConfig = {'FontName', 'Microsoft YaHei', 'FontSize', 14};

fig = figure(figureConfig{:});

% 子图1：各年份工业负荷曲线预测
subplot(1, 2, 1);
hold on;
colors = lines(3); % 获取不同颜色

p1 = plot(t_axis, results.y2024.load / Unit_YiKW_to_KW, 'Color', colors(1,:), 'LineWidth', 2, 'Marker', 'o', 'MarkerIndices', 1:3:24);
p2 = plot(t_axis, results.y2025.load / Unit_YiKW_to_KW, 'Color', colors(2,:), 'LineWidth', 2, 'Marker', 's', 'MarkerIndices', 1:3:24);
p3 = plot(t_axis, results.y2030.load / Unit_YiKW_to_KW, 'Color', colors(3,:), 'LineWidth', 2, 'Marker', '^', 'MarkerIndices', 1:3:24);

title('工业大负荷时序预测曲线 (2024-2030)', 'FontWeight', 'bold', 'FontSize', 16);
xlabel('时间 (小时)', fontConfig{:});
ylabel('工业负荷 (亿千瓦)', fontConfig{:});
legend([p1, p2, p3], {'2024年 (1.145亿基准)', '2025年 (1.302亿基准)', '2030年 (1.6亿基准)'}, 'Location', 'northwest', fontConfig{:});
grid on;
ax = gca; ax.FontSize = 12;
xlim([0 23]);
xticks(0:2:23);

% 子图2：调节能力预测
subplot(1, 2, 2);
hold on;

r1 = plot(t_axis, results.y2024.reg / Unit_WanKW_to_KW, 'Color', colors(1,:), 'LineWidth', 2, 'LineStyle', '--');
r2 = plot(t_axis, results.y2025.reg / Unit_WanKW_to_KW, 'Color', colors(2,:), 'LineWidth', 2, 'LineStyle', '--');
r3 = plot(t_axis, results.y2030.reg / Unit_WanKW_to_KW, 'Color', colors(3,:), 'LineWidth', 2, 'LineStyle', '--');

title('工业负荷可调节潜力预测 (按5%测算)', 'FontWeight', 'bold', 'FontSize', 16);
xlabel('时间 (小时)', fontConfig{:});
ylabel('调节能力 (万千瓦)', fontConfig{:});
legend([r1, r2, r3], {'2024年调节潜力', '2025年调节潜力', '2030年调节潜力'}, 'Location', 'northwest', fontConfig{:});
grid on;
ax = gca; ax.FontSize = 12;
xlim([0 23]);
xticks(0:2:23);

% 保存图片
print(fig, 'Industrial_Prediction_2024_2030.png', '-dpng', '-r300');