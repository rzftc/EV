function P_base_sequence = EVbaseP_Average(C_EV, eta, E_tar, E_in, ...
                                           t_dep, t_in, dt, varargin)
% EVbaseP_Average 计算“平均/线性充电”策略的基线功率
% 输入参数与 EVbaseP_ChargeUntilFull 保持一致以便替换

    % 提取额外参数 (varargin) 以兼容调用格式
    % 这里的 p_on (额定功率) 仅作为上限约束
    if length(varargin) < 5
        error('输入参数不足');
    end
    p_on_max = varargin{2};     % 额定最大功率
    num_time_points = varargin{4}; 
    time_axis = varargin{5};    % 时间轴

    P_base_sequence = zeros(num_time_points, 1);
    
    % 计算需要的平均功率
    duration_h = t_dep - t_in;
    if duration_h <= 0
        P_avg = 0;
    else
        energy_needed = E_tar - E_in;
        P_avg = energy_needed / (eta * duration_h);
    end
    
    % 功率必须在 [0, p_on_max] 之间
    % 如果需求超过额定功率，则只能全速充（此时无上调潜力，但这是物理限制）
    P_base_applied = min(max(P_avg, 0), p_on_max);

    % 生成序列
    for t_idx = 1:num_time_points
        current_t = time_axis(t_idx);
        % 判断是否在线
        if current_t >= t_in && current_t < t_dep
            P_base_sequence(t_idx) = P_base_applied;
        else
            P_base_sequence(t_idx) = 0;
        end
    end
end