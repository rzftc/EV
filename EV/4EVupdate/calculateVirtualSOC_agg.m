function [P_agg, S_agg_next] = calculateVirtualSOC_agg(EVs, dt)
    % calculateVirtualSOC_agg: 计算聚合体的总功率和下一时刻的聚合期望SOC
    % 
    % 修改说明：
    % 为了保证 S_agg 与 lambda (控制信号) 的物理含义一致，
    % 此处筛选车辆的逻辑必须与 aggregateEVs.m 中 normal_EVs 的筛选逻辑完全同步。
    % 即：剔除 LockON (刚性负荷) 和 LockOFF (不在线/已满)，仅计算“可调节资源”的状态。

    %% 1. 车辆筛选 (关键修改)
    % 逻辑：active_mask = NOT LockON AND NOT LockOFF
    % 在 updateLockState 中，主要状态只有 ON, LockON, LockOFF。
    % 因此这实际上等同于只选择 state == 'ON' 的车辆。
    
    lockON_mask  = strcmp({EVs.state}, 'LockON');
    lockOFF_mask = strcmp({EVs.state}, 'LockOFF');
    
    % 仅保留可参与调节的“正常状态”车辆
    active_mask = ~lockON_mask & ~lockOFF_mask; 
    active_EVs = EVs(active_mask);
    
    %% 2. 边界条件处理 (防止除以零导致 NaN)
    if isempty(active_EVs)
        P_agg = 0;
        % 当没有可控车辆时，S_agg 定义为 0 (中性) 或保持上一时刻值
        % 这里返回 0 以避免 NaN 污染绘图
        S_agg_next = 0; 
        return;
    end
    
    %% 3. 聚合状态计算
    P_agg = 0;
    P_req = 0;
    E_exp = 0;
    E_actual = 0;
    C = 0;
    
    % 预计算时间步长 (小时)
    dt_h = dt / 60;
    
    % 遍历有效车辆进行累加
    for i = 1:length(active_EVs) 
        % 累加实时功率
        P_agg = P_agg + active_EVs(i).P_current;
        
        % 计算该车辆的期望充电功率 P_req (理想轨迹斜率)
        time_window = (active_EVs(i).t_dep - active_EVs(i).t_in) / 60;
        if time_window > 0
            P_req = (active_EVs(i).E_tar - active_EVs(i).E_ini) / (active_EVs(i).eta * time_window);
        else
            P_req = 0;
        end
        
        % 累加期望能量 E_exp (理想轨迹)
        % 注意：这里使用 active_EVs(i).E_exp 是累积量，加上当前步的变化
        E_exp = E_exp + active_EVs(i).E_exp + active_EVs(i).eta * P_req * dt_h;
        
        % 累加实际能量 E_actual
        E_actual = E_actual + active_EVs(i).E_actual + active_EVs(i).eta * active_EVs(i).P_current * dt_h;
        
        % 累加容量 (作为加权分母)
        C = C + active_EVs(i).C;
    end
    
    %% 4. 计算 S_agg_next
    % 获取调节系数 r (假设同构或取第一个作为代表，若 r 为异构需改为加权平均)
    r_representative = active_EVs(1).r; 
    
    if C > 0 && r_representative ~= 0
        % 聚合公式：S_agg = - (Total_E_actual - Total_E_exp) / (Total_C * r)
        S_agg_next = -(E_actual - E_exp) / (C * r_representative);
    else
        S_agg_next = 0;
    end
end