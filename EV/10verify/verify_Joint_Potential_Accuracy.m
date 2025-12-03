%% verify_Joint_Potential_Accuracy.m
% 功能：
% 1. 读取 EV 和 AC 的 5分钟仿真结果文件。
% 2. 绘制 EV、AC、联合系统 的 聚合模型 vs 单体累加 对比图 (8张图，独立保存，无标题)。
% 3. 计算并输出 EV与AC联合系统 的上、下调节潜力总量偏差。

clc; clear; close all;

set(0, 'DefaultAxesFontName', 'Microsoft YaHei');
set(0, 'DefaultTextFontName', 'Microsoft YaHei');

%% 1. 文件路径定义
ev_file = 'DER1.mat';
ac_file = 'DER2.mat';

fprintf('==================================================\n');
fprintf('      源荷聚合调节潜力 联合验证程序      \n');
fprintf('==================================================\n');

%% 2. 加载 EV 数据
fprintf('正在加载 DER1 数据: %s ...\n', ev_file);
if ~exist(ev_file, 'file')
    error('未找到 EV 数据文件: %s', ev_file);
end
data_ev = load(ev_file);
if ~isfield(data_ev, 'results')
    error('EV 文件中未找到 "results" 结构体。');
end
ev_res = data_ev.results;

% 提取 EV 数据
EV_Up_Model   = ev_res.EV_Up(:);                   % 聚合模型值
EV_Down_Model = ev_res.EV_Down(:);
EV_Up_True    = ev_res.EV_Up_Individual_Sum(:);    % 单体累加真值
EV_Down_True  = ev_res.EV_Down_Individual_Sum(:);

% 提取时间轴
if isfield(ev_res, 'time_points_absolute')
    time_axis = ev_res.time_points_absolute(:);
else
    dt_short = 5; 
    time_axis = (0:length(EV_Up_Model)-1)' * dt_short / 60;
end

%% 3. 加载 AC 数据
fprintf('正在加载 DER2 数据: %s ...\n', ac_file);
if ~exist(ac_file, 'file')
    error('未找到 AC 数据文件: %s', ac_file);
end
data_ac = load(ac_file);
if ~isfield(data_ac, 'results')
    error('AC 文件中未找到 "results" 结构体。');
end
ac_res = data_ac.results;

% 提取 AC 数据
if isfield(ac_res, 'Agg_Model_Potential_Up_History')
    AC_Up_Model = ac_res.Agg_Model_Potential_Up_History(:);
    AC_Down_Model = ac_res.Agg_Model_Potential_Down_History(:);
else
    error('AC 结果中缺少聚合模型数据 (Agg_Model_Potential_Up_History)。');
end

if isfield(ac_res, 'Agg_P_Potential_Up_History')
    AC_Up_True = ac_res.Agg_P_Potential_Up_History(:);
    AC_Down_True = ac_res.Agg_P_Potential_Down_History(:);
else
    error('AC 结果中缺少单体累加数据 (Agg_P_Potential_Up_History)。');
end

%% 4. 数据对齐
len_min = min([length(EV_Up_Model), length(AC_Up_Model)]);

if length(EV_Up_Model) ~= len_min
    fprintf('提示: 数据长度不一致，截取前 %d 个时间步。\n', len_min);
end

% 截取
EV_Up_Model   = EV_Up_Model(1:len_min);
EV_Down_Model = EV_Down_Model(1:len_min);
EV_Up_True    = EV_Up_True(1:len_min);
EV_Down_True  = EV_Down_True(1:len_min);

AC_Up_Model   = AC_Up_Model(1:len_min);
AC_Down_Model = AC_Down_Model(1:len_min);
AC_Up_True    = AC_Up_True(1:len_min);
AC_Down_True  = AC_Down_True(1:len_min);

time_axis = time_axis(1:len_min);

%% 5. 计算联合系统数据
% 联合模型 = EV模型 + AC模型
Total_Up_Model   = EV_Up_Model + AC_Up_Model;
Total_Down_Model = EV_Down_Model + AC_Down_Model;

% 联合真值 = EV真值 + AC真值
Total_Up_True    = EV_Up_True + AC_Up_True;
Total_Down_True  = EV_Down_True + AC_Down_True;

% 计算比值 (真值 / 聚合值)
epsilon = 1e-3;

Ratio_Up = ones(size(Total_Up_Model));
valid_idx_up = abs(Total_Up_Model) > epsilon;
Ratio_Up(valid_idx_up) = Total_Up_True(valid_idx_up) ./ Total_Up_Model(valid_idx_up);

Ratio_Down = ones(size(Total_Down_Model));
valid_idx_down = abs(Total_Down_Model) > epsilon;
Ratio_Down(valid_idx_down) = Total_Down_True(valid_idx_down) ./ Total_Down_Model(valid_idx_down);

%% 6. 独立绘图与保存
dpi_value = 300; % 高DPI设置

% ================= [1. EV 上调对比] =================
fig1 = figure('Color', 'w', 'Position', [100, 100, 800, 500]);
plot(time_axis, EV_Up_Model, 'r-', 'LineWidth', 1.5, 'DisplayName', '聚合模型');
hold on;
plot(time_axis, EV_Up_True, 'k--', 'LineWidth', 1.5, 'DisplayName', '单体累加');
xlabel('时间 (小时)', 'FontSize', 14); ylabel('功率 (kW)', 'FontSize', 14);
legend('Location', 'best', 'FontSize', 12); grid on; set(gca, 'FontSize', 12);
% 无标题
save_name1 = '电动汽车上调潜力对比.png';
print(fig1, save_name1, '-dpng', sprintf('-r%d', dpi_value));
fprintf('已保存: %s\n', save_name1);

% ================= [2. EV 下调对比] =================
fig2 = figure('Color', 'w', 'Position', [150, 150, 800, 500]);
plot(time_axis, EV_Down_Model, 'b-', 'LineWidth', 1.5, 'DisplayName', '聚合模型');
hold on;
plot(time_axis, EV_Down_True, 'k--', 'LineWidth', 1.5, 'DisplayName', '单体累加');
xlabel('时间 (小时)', 'FontSize', 14); ylabel('功率 (kW)', 'FontSize', 14);
legend('Location', 'best', 'FontSize', 12); grid on; set(gca, 'FontSize', 12);
% 无标题
save_name2 = '电动汽车下调潜力对比.png';
print(fig2, save_name2, '-dpng', sprintf('-r%d', dpi_value));
fprintf('已保存: %s\n', save_name2);

% ================= [3. AC 上调对比] =================
fig3 = figure('Color', 'w', 'Position', [200, 200, 800, 500]);
plot(time_axis, AC_Up_Model, 'r-', 'LineWidth', 1.5, 'DisplayName', '聚合模型');
hold on;
plot(time_axis, AC_Up_True, 'k--', 'LineWidth', 1.5, 'DisplayName', '单体累加');
xlabel('时间 (小时)', 'FontSize', 14); ylabel('功率 (kW)', 'FontSize', 14);
legend('Location', 'best', 'FontSize', 12); grid on; set(gca, 'FontSize', 12);
% 无标题
save_name3 = '空调上调潜力对比.png';
print(fig3, save_name3, '-dpng', sprintf('-r%d', dpi_value));
fprintf('已保存: %s\n', save_name3);

% ================= [4. AC 下调对比] =================
fig4 = figure('Color', 'w', 'Position', [250, 250, 800, 500]);
plot(time_axis, AC_Down_Model, 'b-', 'LineWidth', 1.5, 'DisplayName', '聚合模型');
hold on;
plot(time_axis, AC_Down_True, 'k--', 'LineWidth', 1.5, 'DisplayName', '单体累加');
xlabel('时间 (小时)', 'FontSize', 14); ylabel('功率 (kW)', 'FontSize', 14);
legend('Location', 'best', 'FontSize', 12); grid on; set(gca, 'FontSize', 12);
% 无标题
save_name4 = '空调下调潜力对比.png';
print(fig4, save_name4, '-dpng', sprintf('-r%d', dpi_value));
fprintf('已保存: %s\n', save_name4);

% ================= [5. 联合系统 上调对比] =================
fig5 = figure('Color', 'w', 'Position', [300, 300, 800, 500]);
plot(time_axis, Total_Up_Model, 'r-', 'LineWidth', 1.5, 'DisplayName', '聚合模型');
hold on;
plot(time_axis, Total_Up_True, 'k--', 'LineWidth', 1.5, 'DisplayName', '单体累加');
xlabel('时间 (小时)', 'FontSize', 14); ylabel('功率 (kW)', 'FontSize', 14);
legend('Location', 'best', 'FontSize', 12); grid on; set(gca, 'FontSize', 12);
% 无标题
save_name5 = '联合系统上调潜力对比.png';
print(fig5, save_name5, '-dpng', sprintf('-r%d', dpi_value));
fprintf('已保存: %s\n', save_name5);

% ================= [6. 联合系统 下调对比] =================
fig6 = figure('Color', 'w', 'Position', [350, 350, 800, 500]);
plot(time_axis, Total_Down_Model, 'b-', 'LineWidth', 1.5, 'DisplayName', '聚合模型');
hold on;
plot(time_axis, Total_Down_True, 'k--', 'LineWidth', 1.5, 'DisplayName', '单体累加');
xlabel('时间 (小时)', 'FontSize', 14); ylabel('功率 (kW)', 'FontSize', 14);
legend('Location', 'best', 'FontSize', 12); grid on; set(gca, 'FontSize', 12);
% 无标题
save_name6 = '联合系统下调潜力对比.png';
print(fig6, save_name6, '-dpng', sprintf('-r%d', dpi_value));
fprintf('已保存: %s\n', save_name6);

% ================= [7. 联合系统 上调比值] =================
fig7 = figure('Color', 'w', 'Position', [400, 400, 800, 500]);
plot(time_axis, Ratio_Up, 'm-', 'LineWidth', 1.5);
hold on;
yline(1.0, 'k--', 'LineWidth', 1.5); % 参考线
xlabel('时间 (小时)', 'FontSize', 14); ylabel('比值 (真值/聚合值)', 'FontSize', 14);
ylim([0.5, 1.5]); grid on; set(gca, 'FontSize', 12);
% 无标题
save_name7 = '联合系统上调潜力真值与聚合值之比.png';
print(fig7, save_name7, '-dpng', sprintf('-r%d', dpi_value));
fprintf('已保存: %s\n', save_name7);

% ================= [8. 联合系统 下调比值] =================
fig8 = figure('Color', 'w', 'Position', [450, 450, 800, 500]);
plot(time_axis, Ratio_Down, 'c-', 'LineWidth', 1.5);
hold on;
yline(1.0, 'k--', 'LineWidth', 1.5); % 参考线
xlabel('时间 (小时)', 'FontSize', 14); ylabel('比值 (真值/聚合值)', 'FontSize', 14);
ylim([0.5, 1.5]); grid on; set(gca, 'FontSize', 12);
% 无标题
save_name8 = '联合系统下调潜力真值与聚合值之比.png';
print(fig8, save_name8, '-dpng', sprintf('-r%d', dpi_value));
fprintf('已保存: %s\n', save_name8);


%% 7. 总量偏差计算
fprintf('\n========== 可控资源可调潜力预测验证 ==========\n');

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
fprintf('  - 可控资源与电网时空互济规划态潜力预测平均精确度: %.2f%%\n',      100-Err_Up);

fprintf('【联合下调潜力】\n');
fprintf('  - 聚合模型总量:                                %.2f kW·step\n', Sum_Model_Down);
fprintf('  - 单体累加总量:                                %.2f kW·step\n', Sum_True_Down);
fprintf('  - 可控资源与电网时空互济规划态潜力预测平均精确度: %.2f%%\n',      100-Err_Down);

fprintf('==================================================\n');