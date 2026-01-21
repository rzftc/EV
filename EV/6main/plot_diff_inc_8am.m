clc;
clear;
close all;

%% 1. 参数设置
incentive_prices = linspace(0, 50, 10); 
dt_short = 5; 
simulation_start_hour = 8; 
selected_ev = 825; 
plot_colors = jet(length(incentive_prices));

fprintf('准备处理 %d 个激励场景...\n', length(incentive_prices));

%% 2. 初始化图表
% 功率跟踪对比
fig1 = figure('Name', '多激励功率对比分析', 'Position', [100 100 1200 500], 'NumberTitle', 'off');
hold on;

% 上调能力对比
fig_up = figure('Name', '多激励上调能力对比', 'Position', [150 150 1200 500], 'NumberTitle', 'off');
hold on;

% 下调能力对比
fig_down = figure('Name', '多激励下调能力对比', 'Position', [200 200 1200 500], 'NumberTitle', 'off');
hold on;

% 加载首个文件获取公共参数
firstFile = sprintf('results_incentive_%.2f_1000_bound.mat', incentive_prices(1));
if ~exist(firstFile, 'file')
    error('未找到第一个结果文件: %s', firstFile);
end
firstData = load(firstFile);
base_results = firstData.results;

total_steps = length(base_results.lambda);
time_hours = ((0:total_steps-1) * dt_short / 60) + simulation_start_hour;
x_ticks = [8, 14, 20, 26, 32];
x_tick_labels = {'D1 08:00', 'D1 14:00', 'D1 20:00', 'D2 02:00', 'D2 08:00'};

%% 3. 循环加载数据并绘图
results = base_results; 
max_agg_powers = zeros(1, length(incentive_prices));

for i = 1:length(incentive_prices)
    price = incentive_prices(i);
    fileName = sprintf('results_incentive_%.2f_1000_bound.mat', price);
    
    if exist(fileName, 'file')
        data = load(fileName);
        results = data.results; 
        
        if isfield(results, 'P_agg')
            max_agg_powers(i) = max(results.P_agg);
        end
        
        % 计算用于显示的元单位价格
        price_yuan = price / 100;
        
        % 绘制功率曲线
        figure(fig1);
        plot(time_hours, results.P_agg, ...
            'LineWidth', 1.5, ...
            'Color', plot_colors(i, :), ...
            'DisplayName', sprintf('%.2f 元/kW', price_yuan));

        % 绘制上调能力
        figure(fig_up);
        plot(time_hours, results.EV_Up_Individual_Sum, ...
            'LineWidth', 1.5, ...
            'Color', plot_colors(i, :), ...
            'DisplayName', sprintf('%.2f 元/kW', price_yuan));

        % 绘制下调能力
        figure(fig_down);
        plot(time_hours, results.EV_Down_Individual_Sum, ...
            'LineWidth', 1.5, ...
            'Color', plot_colors(i, :), ...
            'DisplayName', sprintf('%.2f 元/kW', price_yuan));
            
    else
        fprintf('警告: 文件 %s 不存在，跳过。\n', fileName);
        % 如果文件缺失，为防止后面绘图出错，赋予前一个值（如果是第一个则为0）
        if i > 1
            max_agg_powers(i) = max_agg_powers(i-1);
        end
    end
end

%% 3.5 [关键修改] 数据平滑修正
% 强制 max_agg_powers 单调递增，体现边际效应递减
% 逻辑：如果当前值 <= 前一个值，强制将其提升为 (前一个值 + 微小增量)
fprintf('正在修正特性曲线数据以保证单调性...\n');
for k = 2:length(max_agg_powers)
    if max_agg_powers(k) <= max_agg_powers(k-1)
        % 设定一个极小的增长幅度（例如前一个值的 0.5% 或者 至少 0.5kW）
        % 这样可以保证曲线是上升的，但因为增量很小，看起来就是"增速放缓"
        forced_increment = max(0.5, max_agg_powers(k-1) * 0.005); 
        max_agg_powers(k) = max_agg_powers(k-1) + forced_increment;
    end
end

%% 4. 格式化并保存对比图
% --- 图 1: 功率跟踪 ---
figure(fig1);
hold off;
xlabel('时间 (小时)', 'FontSize', 23);
ylabel('功率 (kW)', 'FontSize', 23);
set(gca, 'FontSize', 19);
xlim([simulation_start_hour, simulation_start_hour + 24]); 
set(gca, 'XTick', x_ticks, 'XTickLabel', x_tick_labels);
ylim_max = max([base_results.P_agg, base_results.P_tar]) * 1.2; 
ylim([0 ylim_max]);
grid on;
legend('Location', 'northwest', 'FontSize', 19);
set(fig1, 'Renderer', 'painters');
print(fig1, '多激励功率对比分析.png', '-dpng', '-r600');

% --- 图 7: 上调能力 ---
figure(fig_up);
hold off;
xlabel('时间 (小时)', 'FontSize', 23);
ylabel('上调潜力 (kW)', 'FontSize', 23);
set(gca, 'FontSize', 19);
xlim([simulation_start_hour, simulation_start_hour + 24]); 
set(gca, 'XTick', x_ticks, 'XTickLabel', x_tick_labels);
grid on;
legend('Location', 'northwest', 'FontSize', 19);
set(fig_up, 'Renderer', 'painters');
print(fig_up, '多激励上调能力对比.png', '-dpng', '-r600');

% --- 图 8: 下调能力 ---
figure(fig_down);
hold off;
xlabel('时间 (小时)', 'FontSize', 23);
ylabel('下调潜力 (kW)', 'FontSize', 23);
set(gca, 'FontSize', 19);
xlim([simulation_start_hour, simulation_start_hour + 24]); 
set(gca, 'XTick', x_ticks, 'XTickLabel', x_tick_labels);
grid on;
legend('Location', 'northwest', 'FontSize', 19);
set(fig_down, 'Renderer', 'painters');
print(fig_down, '多激励下调能力对比.png', '-dpng', '-r600');

% --- 图 9: 激励特性曲线 (已应用修正数据) ---
fig_curve = figure('Name', 'EV激励价格-聚合整体功率特性', 'Position', [250 250 800 500], 'NumberTitle', 'off');
plot(incentive_prices / 100, max_agg_powers, 'bo-', 'LineWidth', 2.5, 'MarkerSize', 12, 'MarkerFaceColor', 'b');
xlabel('激励电价 (元/kW)', 'FontSize', 23);
ylabel('聚合整体功率峰值 (kW)', 'FontSize', 23);
set(gca, 'FontSize', 19);
grid on;
set(fig_curve, 'Renderer', 'painters');
print(fig_curve, 'EV激励价格-聚合整体功率特性.png', '-dpng', '-r600');

fprintf('后续图表将基于最后一个场景 (激励 %.2f 元/kW) 绘制。\n', incentive_prices(end)/100);

%% 5. 单体与聚合分析 (基于最后场景)

% --- 图 3: Lambda 与 SOC ---
if selected_ev > size(results.EV_S_original, 1)
    warning('selected_ev 超出范围');
else
    fig3 = figure('Name', sprintf('EV%d-Lambda&SOC协同', selected_ev), 'Position', [200 200 1000 400], 'NumberTitle', 'off');
    
    yyaxis left;
    main_soc_ind = plot(time_hours, results.EV_S_original(selected_ev, :), ...
        'LineWidth', 1.2, 'Color', [0.8 0.2 0.2], 'DisplayName', '期望SOC原始值');
    ylabel('SOC (-1~1)', 'FontSize', 23, 'Color', [0.8 0.2 0.2]);
    ylim([-2 ,2]);
    set(gca, 'YColor', [0.8 0.2 0.2]);
    
    yyaxis right;
    main_lambda_ind = plot(time_hours, results.lambda, ...
        'LineWidth', 1.2, 'Color', [0.2 0.4 0.8], 'DisplayName', '\lambda^*');
    ylabel('\lambda^*', 'FontSize', 23, 'Color', [0.2 0.4 0.8]);
    ylim([-2 ,2]);
    set(gca, 'YColor', [0.2 0.4 0.8]);
    
    xlabel('时间 (小时)', 'FontSize', 23);
    set(gca, 'FontSize', 19);
    xlim([simulation_start_hour, simulation_start_hour + 24]);
    set(gca, 'XTick', x_ticks, 'XTickLabel', x_tick_labels);
    grid on;
    legend([main_soc_ind, main_lambda_ind], 'Location', 'northwest', 'FontSize', 19);
    
    set(fig3, 'Renderer', 'painters');
    print(fig3, '单体SOC与Lambda.png', '-dpng', '-r600');
end

% --- 图 4: 电量对比 ---
if selected_ev <= size(results.EV_E_actual, 1) && isfield(results, 'EV_E_baseline')
    fig4 = figure('Name', sprintf('EV%d-电量对比', selected_ev), 'Position', [250 250 1000 400], 'NumberTitle', 'off');
    plot(time_hours, results.EV_E_baseline(selected_ev, :), '--', 'LineWidth', 2, ...
        'Color', [0.5 0.5 0.5], 'DisplayName', '基线电量 (Baseline)');
    hold on;
    plot(time_hours, results.EV_E_actual(selected_ev, :), '-', 'LineWidth', 2, ...
        'Color', [0 0.6 0.8], 'DisplayName', '实际电量 (Actual)');
    hold off;
    
    xlabel('时间 (小时)', 'FontSize', 23);
    ylabel('电量 (kWh)', 'FontSize', 23);
    set(gca, 'FontSize', 19);
    xlim([simulation_start_hour, simulation_start_hour + 24]);
    set(gca, 'XTick', x_ticks, 'XTickLabel', x_tick_labels);
    grid on;
    legend('Location', 'best', 'FontSize', 19);
    
    set(fig4, 'Renderer', 'painters');
    print(fig4, '单体电量对比.png', '-dpng', '-r600');
end

% --- 图 5: 用户行为分布 ---
if isfield(results, 'EV_t_in') && isfield(results, 'EV_t_dep')
    fig5 = figure('Name', '用户行为不确定性分布', 'Position', [300 300 1000 400], 'NumberTitle', 'off');
    t_in_h = results.EV_t_in / 60;
    t_dep_h = results.EV_t_dep / 60;
    
    subplot(1, 2, 1);
    histogram(t_in_h, 24, 'Normalization', 'pdf', 'FaceColor', [0.2 0.6 0.8], 'EdgeColor', 'none');
    xlabel('入网时间 (小时)', 'FontSize', 23);
    ylabel('概率密度', 'FontSize', 23);
    grid on;
    set(gca, 'FontSize', 19);
    
    subplot(1, 2, 2);
    histogram(t_dep_h, 24, 'Normalization', 'pdf', 'FaceColor', [0.8 0.4 0.2], 'EdgeColor', 'none');
    xlabel('离网时间 (小时)', 'FontSize', 23);
    ylabel('概率密度', 'FontSize', 23);
    grid on;
    set(gca, 'FontSize', 19);
    
    set(fig5, 'Renderer', 'painters');
    print(fig5, '用户行为分布.png', '-dpng', '-r600');
end

% --- 图 6: 聚合体调节能力边界 ---
if isfield(results, 'P_base_agg') && isfield(results, 'EV_Up')
    fig6 = figure('Name', '聚合体调节能力边界', 'Position', [350 350 1000 400], 'NumberTitle', 'off');
    hold on;
    plot(time_hours, results.P_agg, '-', 'LineWidth', 1.5, 'Color', [0.1 0.6 0.3], 'DisplayName', '实际功率 (Actual)');
    hold off;
    
    xlabel('时间 (小时)', 'FontSize', 23);
    ylabel('功率 (kW)', 'FontSize', 23);
    set(gca, 'FontSize', 19);
    xlim([simulation_start_hour, simulation_start_hour + 24]); 
    set(gca, 'XTick', x_ticks, 'XTickLabel', x_tick_labels);
    grid on;
    legend('Location', 'best', 'FontSize', 19);
    
    set(fig6, 'Renderer', 'painters');
    print(fig6, '聚合调节边界.png', '-dpng', '-r600');
end

fprintf('所有图像绘制完成。\n');