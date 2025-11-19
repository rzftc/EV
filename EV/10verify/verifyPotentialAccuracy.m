function verifyPotentialAccuracy(results, dt_short)
    % verifyPotentialAccuracy 验证聚合模型的准确性 (基于总量偏差)
    %
    % 输入:
    %   results  - 包含 EV_Up, EV_Down, EV_Up_Individual_Sum, EV_Down_Individual_Sum 的结构体
    %   dt_short - 时间步长 (用于绘图时间轴)
    
    fprintf('\n==================================================\n');
    fprintf('      EV 聚合调节潜力 准确性验证报告 (总量版)      \n');
    fprintf('==================================================\n');

    %% 1. 数据提取
    % 模型预测值 (Aggregated Model)
    Up_Model   = results.EV_Up;
    Down_Model = results.EV_Down;
    
    % 真值 (Individual Sum)
    Up_True   = results.EV_Up_Individual_Sum;
    Down_True = results.EV_Down_Individual_Sum;
    
    % 时间轴
    if nargin < 2
        dt_short = 5; 
    end
    time_axis = (0:length(Up_Model)-1) * dt_short / 60;
    
    %% 2. 总量误差计算 (Total Capacity Error)
    % 计算全时段的总潜力 (相当于对功率积分，忽略dt常数不影响百分比)
    
    % --- 上调潜力 (Up) ---
    Total_Up_Model = sum(Up_Model);
    Total_Up_True  = sum(Up_True);
    
    % 防止分母为0
    if abs(Total_Up_True) > 1e-3
        % 计算总量的相对误差
        Err_Up_Total_Pct = abs(Total_Up_Model - Total_Up_True) / abs(Total_Up_True) * 100;
    else
        Err_Up_Total_Pct = 0;
        fprintf('提示: 上调潜力总量接近 0，跳过误差计算。\n');
    end
    
    % --- 下调潜力 (Down) ---
    Total_Down_Model = sum(Down_Model);
    Total_Down_True  = sum(Down_True);
    
    if abs(Total_Down_True) > 1e-3
        Err_Down_Total_Pct = abs(Total_Down_Model - Total_Down_True) / abs(Total_Down_True) * 100;
    else
        Err_Down_Total_Pct = 0;
        fprintf('提示: 下调潜力总量接近 0，跳过误差计算。\n');
    end
    
    %% 3. 结果输出与判定
    fprintf('\n【上调潜力 (Up-Regulation)】\n');
    fprintf('  - 模型总量:       %.2f kW·step\n', Total_Up_Model);
    fprintf('  - 真值总量:       %.2f kW·step\n', Total_Up_True);
    fprintf('  - 总量偏差:       %.2f%%\n', Err_Up_Total_Pct);
    
    if Err_Up_Total_Pct < 10
        fprintf('  >>> 结果: [合格] (偏差 < 10%%)\n');
    else
        fprintf('  >>> 结果: [警告] (偏差 >= 10%%)\n');
    end
    
    fprintf('\n【下调潜力 (Down-Regulation)】\n');
    fprintf('  - 模型总量:       %.2f kW·step\n', Total_Down_Model);
    fprintf('  - 真值总量:       %.2f kW·step\n', Total_Down_True);
    fprintf('  - 总量偏差:       %.2f%%\n', Err_Down_Total_Pct);
    
    if Err_Down_Total_Pct < 10
        fprintf('  >>> 结果: [合格] (偏差 < 10%%)\n');
    else
        fprintf('  >>> 结果: [警告] (偏差 >= 10%%)\n');
    end
    
    %% 4. 可视化对比
    % 绘制累积潜力曲线 (Cumulative Potential)，直观展示总量差异
    figure('Name', '聚合模型总量累积对比', 'Position', [200, 200, 1000, 400]);
    
    subplot(1, 2, 1);
    plot(time_axis, cumsum(Up_Model), 'r-', 'LineWidth', 1.5, 'DisplayName', '模型累积值');
    hold on;
    plot(time_axis, cumsum(Up_True), 'k--', 'LineWidth', 1.5, 'DisplayName', '真值累积值');
    title(sprintf('上调潜力累积对比 (误差: %.1f%%)', Err_Up_Total_Pct));
    xlabel('时间 (小时)'); ylabel('累积容量'); legend('Location', 'best'); grid on;
    
    subplot(1, 2, 2);
    plot(time_axis, cumsum(Down_Model), 'b-', 'LineWidth', 1.5, 'DisplayName', '模型累积值');
    hold on;
    plot(time_axis, cumsum(Down_True), 'k--', 'LineWidth', 1.5, 'DisplayName', '真值累积值');
    title(sprintf('下调潜力累积对比 (误差: %.1f%%)', Err_Down_Total_Pct));
    xlabel('时间 (小时)'); ylabel('累积容量'); legend('Location', 'best'); grid on;
    
    fprintf('\n==================================================\n');
end