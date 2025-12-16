%% 分布式光伏时空互济潜力预测主程序
% 基于技术报告第3.2.4节 分布式光伏模型 (公式3-63 ~ 3-74)
% 数据来源: 附件2 (2024年全省风光负荷曲线) & 用户提供的规划数据
% 修改时间: 2025-12-08

clear; clc; close all;

%% 1. 基础参数设置 (万千瓦 -> kW)
Unit_WanKW_to_KW = 10^4;

% 规划装机容量 (用户输入)
Cap_2025_Oct = 6343 * Unit_WanKW_to_KW; % 2025年10月装机: 6343万千瓦
Cap_2030     = 7000 * Unit_WanKW_to_KW; % 2030年预计装机: 7000万千瓦

fprintf('---------- 分布式光伏规模设定 ----------\n');
fprintf('2025年10月装机目标: %.2f 万千瓦\n', Cap_2025_Oct / Unit_WanKW_to_KW);
fprintf('2030年预计装机目标: %.2f 万千瓦\n', Cap_2030 / Unit_WanKW_to_KW);
fprintf('--------------------------------------\n');

%% 2. 读取附件数据 (提取出力特性)
% 文件名
filename = '附件2 2024年全省风光负荷曲线.xlsx';

% 检查文件
if ~isfile(filename)
    % 尝试在上一级目录查找
    if isfile(['../' filename])
        filename = ['../' filename];
    else
        error('错误: 未找到文件 "%s"，请确认路径。', filename);
    end
end

% 读取数据 (保留原始变量名以便索引)
opts = detectImportOptions(filename);
opts.VariableNamingRule = 'preserve';
data_table = readtable(filename, opts);

% 提取关键列
% 假设CSV中包含 '月', '日', '时', '光伏出力率' 列
% 根据附件内容，列名可能为 "月", "日", "光伏出力率"
try
    Months = data_table.('月');
    Days = data_table.('日');
    Hours = data_table.('时');
    PV_Output_Rate = data_table.('光伏出力率'); % 归一化系数 (0~1)
catch
    error('错误: 无法在CSV中找到所需的列名，请检查附件表头。');
end

%% 3. 场景A: 2025年10月全月模拟
% 筛选10月份数据
idx_oct = (Months == 10);
PV_Rate_Oct = PV_Output_Rate(idx_oct);
Time_Oct = 1:length(PV_Rate_Oct); % 简单时间轴

% 计算功率曲线 (P = P_rated * eta)
% 对应技术报告聚合模型思想
P_PV_2025_Oct = Cap_2025_Oct * PV_Rate_Oct;

% 统计特征
Max_P_2025 = max(P_PV_2025_Oct);
Mean_P_2025 = mean(P_PV_2025_Oct);

fprintf('\n[仿真结果] 2025年10月:\n');
fprintf('  最大出力: %.2f 万千瓦\n', Max_P_2025 / Unit_WanKW_to_KW);
fprintf('  平均出力: %.2f 万千瓦\n', Mean_P_2025 / Unit_WanKW_to_KW);

%% 4. 场景B: 2030年典型日模拟 (极限出力场景)
% 选取全年光伏出力率最大的一天作为典型日，评估电网最大消纳压力
[~, idx_max_year] = max(PV_Output_Rate);
target_month = Months(idx_max_year);
target_day = Days(idx_max_year);

% 提取该日24小时数据
idx_typical = (Months == target_month) & (Days == target_day);
PV_Rate_Typical = PV_Output_Rate(idx_typical);
t_axis_day = 0:23;

% 计算2030年典型日功率
P_PV_2030_Typical = Cap_2030 * PV_Rate_Typical;

fprintf('\n[仿真结果] 2030年典型极限日 (%d月%d日):\n', target_month, target_day);
fprintf('  峰值出力: %.2f 万千瓦\n', max(P_PV_2030_Typical) / Unit_WanKW_to_KW);

%% 5. 结果可视化
figure('Color', 'white', 'Position', [100, 100, 1200, 800]);

% 子图1: 2025年10月全月出力曲线
subplot(2, 1, 1);
plot(Time_Oct, P_PV_2025_Oct / Unit_WanKW_to_KW, 'b-', 'LineWidth', 1);
title('2025年10月山东省分布式光伏出力预测 (装机6343万千瓦)', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('时间 (小时)');
ylabel('光伏出力 (万千瓦)');
grid on;
xlim([1, length(Time_Oct)]);
% 添加说明
text(10, max(P_PV_2025_Oct)/Unit_WanKW_to_KW * 0.9, ...
    sprintf('基于2024年同期气象数据外推\n峰值: %.1f 万kW', Max_P_2025/Unit_WanKW_to_KW), ...
    'EdgeColor', 'k', 'BackgroundColor', 'w');

% 子图2: 2030年典型日出力对比
subplot(2, 1, 2);
hold on;
% 绘制2025年同日对比 (假设同气象条件)
P_PV_2025_Typical = Cap_2025_Oct * PV_Rate_Typical;
plot(t_axis_day, P_PV_2025_Typical / Unit_WanKW_to_KW, 'b--', 'LineWidth', 1.5, 'DisplayName', '2025规模 (6343万kW)');
% 绘制2030年预测
area(t_axis_day, P_PV_2030_Typical / Unit_WanKW_to_KW, 'FaceColor', [1 0.8 0.4], 'EdgeColor', 'r', 'LineWidth', 1.5, 'FaceAlpha', 0.5, 'DisplayName', '2030规模 (7000万kW)');

title(sprintf('典型大发日出力对比 (选取自2024年%d月%d日数据)', target_month, target_day), 'FontSize', 12, 'FontWeight', 'bold');
xlabel('时间 (小时)');
ylabel('光伏出力 (万千瓦)');
legend('Location', 'best');
grid on;
xticks(0:2:24);
xlim([0 23]);

% 标记峰值差
peak_2025 = max(P_PV_2025_Typical)/Unit_WanKW_to_KW;
peak_2030 = max(P_PV_2030_Typical)/Unit_WanKW_to_KW;
diff_val = peak_2030 - peak_2025;
text(12, peak_2030, sprintf('增量: +%.1f 万kW', diff_val), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 10, 'Color', 'r');

% 保存图片
print(gcf, 'PV_Prediction_2025_2030.png', '-dpng', '-r300');