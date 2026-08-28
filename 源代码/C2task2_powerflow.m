function [V, delta, P_gen, Q_gen, iter,Y] = C2task2_powerflow()
    % 牛顿-拉夫逊法潮流计算程序
    % 针对课程设计C2的电网结构（已增加节点7为发电机内部节点）
    % 修改：加入了变压器变比模型
    clc;

    %% 1. 系统基准值
    S_base = 1230;               % MVA
    V_base_230 = 230;            % kV (230kV侧基准电压)
    V_base_10_5 = 10.5;          % kV (10.5kV侧基准电压)
    Z_base_230 = V_base_230^2 / S_base;  % Ω
    Z_base_10_5 = V_base_10_5^2 / S_base;% Ω

    %% 2. 节点数据
    n_bus = 7;
    bus_data = [
        3, 1.00, 0,   0,   0;    % 1 平衡
        1, 1.00, 0, 207, 108;    % 2 PQ
        1, 1.00, 0, 315, 171;    % 3 PQ
        1, 1.00, 0, 225, 76.5;    % 4 PQ
        1, 1.00, 0, 360, 175.5;    % 5 PQ
        1, 1.00, 0,   0,   0;    % 6 PQ（高压母线）
        2, 1.05, 0,   0,   0     % 7 PV（发电机内节点）
    ];
    bus_data(:,4:5) = bus_data(:,4:5) / S_base;   % 负荷转标幺

    %% 3. 发电机数据（接在节点7）
    gen_pu = 810*(1-0.07)/S_base;   % 净出力753.3 MW   
    gen_data = [7, gen_pu, 1.05];   % 节点7，PV

    %% 4. 线路数据（包含变压器变比）
    % 单位长度参数→pu/km (基于230kV侧)
    r_pu = 0.054  / Z_base_230;
    x_pu = 0.308  / Z_base_230;
    b_pu = 0.0116e-6 * 2*pi*50 * Z_base_230;   % 2.796e-4 pu/km

    % 线路长度
    line_lengths = [30, 17, 6, 18, 5, 24];   % 对应1-2,2-3,2-4,3-6,4-5,6-5
    
    % 变压器变比计算 (10.5/242 → 标幺变比)
    % 标幺变比 t_pu = (实际变比) × (V_base_to / V_base_from)
    % 对于变压器6→7: from=6(230kV), to=7(10.5kV)
    % 实际变比 = 242/10.5 ≈ 23.0476
    % t_pu = (242/10.5) × (10.5/230) = 242/230 = 1.0522
    t_67 = 242/230;  % 标幺变比

    % 线路参数矩阵 [from, to, r, x, b, t]
    % 普通线路: t=1; 变压器支路: t≠1
    line_data = [
        1, 2, r_pu*30, x_pu*30, b_pu*30/2, 1;
        1, 2, r_pu*30, x_pu*30, b_pu*30/2, 1;
        2, 3, r_pu*17, x_pu*17, b_pu*17/2, 1;
        2, 3, r_pu*17, x_pu*17, b_pu*17/2, 1;
        2, 4, r_pu*6,  x_pu*6,  b_pu*6/2,  1;
        2, 4, r_pu*6,  x_pu*6,  b_pu*6/2,  1;
        3, 6, r_pu*18, x_pu*18, b_pu*18/2, 1;
        3, 6, r_pu*18, x_pu*18, b_pu*18/2, 1;
        4, 5, r_pu*5,  x_pu*5,  b_pu*5/2,  1;
        4, 5, r_pu*5,  x_pu*5,  b_pu*5/2,  1;
        6, 5, r_pu*24, x_pu*24, b_pu*24/2, 1;
        6, 5, r_pu*24, x_pu*24, b_pu*24/2, 1;
        % 变压器支路 6→7
        6, 7, 0, 0.105 * (S_base / 1050), 0, t_67   % X_T = 0.105*1230/1050 ≈ 0.123 pu
    ];

    %% 5. 形成Y矩阵（考虑变比）
    Y = form_y_matrix_with_taps(line_data, n_bus);

    %% 6. 初始化
    V     = bus_data(:,2);
    delta = bus_data(:,3);
    P_load = bus_data(:,4);
    Q_load = bus_data(:,5);
    P_gen  = zeros(n_bus,1);
    Q_gen  = zeros(n_bus,1);
    for k = 1:size(gen_data,1)
        idx = gen_data(k,1);
        P_gen(idx) = gen_data(k,2);
        V(idx)     = gen_data(k,3);
    end
    P_inj = P_gen - P_load;
    Q_inj = Q_gen - Q_load;
    %% 7. 牛顿-拉夫逊迭代（原主循环不变）
    max_iter = 50;  tol = 1e-7;  iter = 0;  converged = false;
    pq_buses = find(bus_data(:,1)==1);
    pv_buses = find(bus_data(:,1)==2);
    slack_bus = find(bus_data(:,1)==3);
    % 无功限值（节点7）
    S_rated = 1050/S_base;
    Q_max = sqrt(S_rated^2 - gen_pu^2);
    Q_min = -Q_max;
    fprintf('发电厂无功限制: Q_min = %.4f pu (%.2f Mvar), Q_max = %.4f pu (%.2f Mvar)\n',...
            Q_min, Q_min*S_base, Q_max, Q_max*S_base);

    % 逻辑索引
    keep_dP = true(n_bus,1); keep_dP(slack_bus) = false;
    keep_dQ = true(n_bus,1); keep_dQ(slack_bus) = false; keep_dQ(pv_buses) = false;

    while iter < max_iter && ~converged
        iter = iter+1;
        [P_calc, Q_calc] = calculate_power(V,delta,Y,n_bus);

        % 检查节点7无功越限
        for i = 1:length(pv_buses)
            bus_idx = pv_buses(i);
            % ------------------------------------------------------------------
            %  节点7无功越限检查 & 类型转换
            for i = 1:length(pv_buses)
              bus_idx = pv_buses(i);
              if Q_calc(bus_idx) > Q_max || Q_calc(bus_idx) < Q_min
               % 1. 决定锁定值
               if Q_calc(bus_idx) > Q_max
                  Q_lock = Q_max;
               else
                  Q_lock = Q_min;
               end
               % 2. 立即锁死注入
               Q_gen(bus_idx) = Q_lock;
               Q_inj(bus_idx) = Q_lock - Q_load(bus_idx);   % 关键！
               % 3. 节点类型转换
               bus_data(bus_idx,1) = 1;           % PV→PQ
               pq_buses = [pq_buses; bus_idx];
               pv_buses(i) = [];
               keep_dQ(bus_idx) = true;
               fprintf('Iteration %d: Bus %d Q越限 (%.3f→%.3f pu) 转为PQ节点\n',...
                     iter, bus_idx, Q_calc(bus_idx), Q_lock);
               % 4. 重新计算功率 & 不平衡量（用锁定值）
               [P_calc, Q_calc] = calculate_power(V,delta,Y,n_bus);
               dP = P_inj - P_calc;
               dQ = Q_inj - Q_calc;
               % 5. 直接构造 mismatch 进入下一步，避免用旧 Q_calc
               dP_red = dP(keep_dP);
               dQ_red = dQ(keep_dQ);
               mismatch = [dP_red; dQ_red];
               if max(abs(mismatch)) < 1e-7
                 converged = true;
               end
               break   % 一次只转一个节点
              end
            end
        end

        dP = P_inj - P_calc;
        dQ = Q_inj - Q_calc;
        dP_red = dP(keep_dP);
        dQ_red = dQ(keep_dQ);
        mismatch = [dP_red; dQ_red];
        if max(abs(mismatch)) < tol, converged = true; break; end

        J = form_jacobian(V,delta,Y,bus_data,n_bus,pq_buses,pv_buses,slack_bus);
        dx = J \ mismatch;
        n_angle = length(dP_red);
        d_angle = dx(1:n_angle);
        dV      = dx(n_angle+1:end);
        angle_idx = find(keep_dP);
        delta(angle_idx) = delta(angle_idx) + d_angle;
        V(pq_buses) = V(pq_buses) + dV;
    end

    %% 8. 结果输出
    if converged
        fprintf('潮流计算在%d次迭代后收敛。\n',iter);
        [P_calc,Q_calc] = calculate_power(V,delta,Y,n_bus);
        calculate_line_flows(V,delta,Y,line_data,S_base);
        display_results(V,delta,P_calc,Q_calc,P_load,Q_load,S_base);
        
        % 正确计算返回的发电功率
        % 计算纯发电功率（节点注入功率 + 负荷功率）
        P_gen_actual = P_calc + P_load;  % pu值
        Q_gen_actual = Q_calc + Q_load;  % pu值
        
        % 转换为实际值(MW, Mvar)返回
        P_gen = P_gen_actual * S_base;
        Q_gen = Q_gen_actual * S_base;
        
    else
        fprintf('潮流计算在%d次迭代后未收敛。\n',iter);
        P_gen = zeros(n_bus,1);
        Q_gen = zeros(n_bus,1);
    end
end
%% 新的Y矩阵形成函数（考虑变比）
function Y = form_y_matrix_with_taps(line_data, n_bus)
    Y = zeros(n_bus, n_bus);
    for i = 1:size(line_data,1)
        from = line_data(i,1);
        to = line_data(i,2);
        r = line_data(i,3);
        x = line_data(i,4);
        b = line_data(i,5);
        t = line_data(i,6); % 变比
        
        z = r + 1j*x;
        y_series = 1/z;
        
        %fprintf("%.4f,%.4f,%.4f\n",r/2,x/2,b*2);
        if abs(t - 1) < 1e-6 % 普通线路
            Y(from, to) = Y(from, to) - y_series;
            Y(to, from) = Y(to, from) - y_series;
            Y(from, from) = Y(from, from) + y_series + 1j*b;
            Y(to, to) = Y(to, to) + y_series + 1j*b;
        else % 变压器支路
            % 标准变压器Π型等值电路
            Y(from, from) = Y(from, from) + y_series/(t^2) + 1j*b/2;
            Y(from, to) = Y(from, to) - y_series/t;
            Y(to, from) = Y(to, from) - y_series/t;
            Y(to, to) = Y(to, to) + y_series + 1j*b/2;
        end
    end
end


function [P_calc,Q_calc] = calculate_power(V,delta,Y,n_bus)
    P_calc = zeros(n_bus,1); Q_calc = zeros(n_bus,1);
    for i = 1:n_bus
        for j = 1:n_bus
            theta = delta(i)-delta(j);
            P_calc(i) = P_calc(i) + V(i)*V(j)*(real(Y(i,j))*cos(theta)+imag(Y(i,j))*sin(theta));
            Q_calc(i) = Q_calc(i) + V(i)*V(j)*(real(Y(i,j))*sin(theta)-imag(Y(i,j))*cos(theta));
        end
    end
end

function J = form_jacobian(V,delta,Y,bus_data,n_bus,pq_buses,pv_buses,slack_bus)
    [P_calc,Q_calc] = calculate_power(V,delta,Y,n_bus);
    n_angle = n_bus-1; n_volt = length(pq_buses);
    J = zeros(n_angle+n_volt,n_angle+n_volt);
    % 计算雅可比矩阵的四个子矩阵: H, N, J, L
    % H = ?P/?δ, N = ?P/?V, J = ?Q/?δ, L = ?Q/?V
    
    % 初始化子矩阵
    H = zeros(n_bus, n_bus);
    N = zeros(n_bus, n_bus);
    J_sub = zeros(n_bus, n_bus);
    L = zeros(n_bus, n_bus);
    
    for i = 1:n_bus
        for j = 1:n_bus
            if i == j
                % 对角线元素
                sum_term = 0;
                for k = 1:n_bus
                    if k ~= i
                        theta_ik = delta(i) - delta(k);
                        sum_term = sum_term + V(k)*(real(Y(i,k))*sin(theta_ik) - imag(Y(i,k))*cos(theta_ik));
                    end
                end
                H(i,i) = -Q_calc(i) - imag(Y(i,i))*V(i)^2;
                N(i,i) = P_calc(i)/V(i) + real(Y(i,i))*V(i);
                J_sub(i,i) = P_calc(i) - real(Y(i,i))*V(i)^2;
                L(i,i) = Q_calc(i)/V(i) - imag(Y(i,i))*V(i);
            else
                % 非对角线元素
                theta_ij = delta(i) - delta(j);
                H(i,j) = V(i)*V(j)*(real(Y(i,j))*sin(theta_ij) - imag(Y(i,j))*cos(theta_ij));
                N(i,j) = V(i)*(real(Y(i,j))*cos(theta_ij) + imag(Y(i,j))*sin(theta_ij));
                J_sub(i,j) = -V(i)*V(j)*(real(Y(i,j))*cos(theta_ij) + imag(Y(i,j))*sin(theta_ij));
                L(i,j) = V(i)*(real(Y(i,j))*sin(theta_ij) - imag(Y(i,j))*cos(theta_ij));
            end
        end
    end
    
    % 移除平衡节点对应的行和列
    all_buses = 1:n_bus;
    non_slack = all_buses(all_buses ~= slack_bus);
    
    % H子矩阵 (?P/?δ)
    H_sub = H(non_slack, non_slack);
    
    % N子矩阵 (?P/?V) - 只保留PQ节点列
    N_sub = N(non_slack, pq_buses);
    
    % J子矩阵 (?Q/?δ) - 只保留PQ节点行
    J_sub2 = J_sub(pq_buses, non_slack);
    
    % L子矩阵 (?Q/?V) - 只保留PQ节点行和列
    L_sub = L(pq_buses, pq_buses);
    
    % 组合成完整的雅可比矩阵
    J = [H_sub, N_sub;
         J_sub2, L_sub];
end

%% 计算线路潮流函数
function calculate_line_flows(V, delta, Y, line_data, S_base)
    fprintf('\n=== 线路潮流分析 ===\n');
    fprintf('线路\t\t有功潮流(MW)\t无功潮流(Mvar)\t有功损耗(MW)\t无功损耗(Mvar)\n');
    
    total_P_loss = 0;
    total_Q_loss = 0;
    
    for i = 1:size(line_data, 1)
        from = line_data(i, 1);
        to = line_data(i, 2);
        r = line_data(i, 3);    % 线路电阻(pu)
        x = line_data(i, 4);    % 线路电抗(pu)
        b = line_data(i, 5);    % 线路电纳(pu)
        
     
        % 计算线路阻抗和导纳
        z = r + 1j*x;
        y_series = 1/z;

        % 计算两端电压相量
        V_from = V(from) * exp(1j*delta(from));
        V_to = V(to) * exp(1j*delta(to));

        % 计算线路电流（从节点流向节点）
        I_from = y_series * (V_from - V_to) + 1j*b * V_from;
        I_to = y_series * (V_to - V_from) + 1j*b * V_to;

        % 计算线路两端功率
        S_from = V_from * conj(I_from) * S_base;
        S_to = V_to * conj(I_to) * S_base;
        % 
        % % 计算线路损耗
        % % I_series = y_series * (V_from - V_to);          % 串联电流 
        % % S_loss = abs(I_series)^2 * z * S_base;          % 真实铜损+铁损
        % % P_loss = real(S_loss);
        % % Q_loss = imag(S_loss);
        P_loss = real(S_from + S_to);
        Q_loss = imag(S_from + S_to);
        % % fprintf('%d→%d  S_from=(%.1f+j%.1f)  S_to=(%.1f+j%.1f)  ΔQ=%.1f\n', ...
        % % from,to, real(S_from),imag(S_from), real(S_to),imag(S_to), Q_loss);

        % 累加总损耗
        total_P_loss = total_P_loss + P_loss;
        total_Q_loss = total_Q_loss + Q_loss;

        % 输出线路潮流结果
        fprintf('%d→%d\t\t%.2f\t\t%.2f\t\t%.2f\t\t%.2f\n', ...
                from, to, real(S_from), imag(S_from), P_loss, Q_loss);
    end
    %% 6→7 变压器损耗明细输出
    % 找到变压器支路（6→7）
    tfIdx = find(line_data(:,1)==6 & line_data(:,2)==7);
    if ~isempty(tfIdx)
        r  = line_data(tfIdx,3);      % 已折算到 1230 MVA 基值
        x  = line_data(tfIdx,4);
        z  = r + 1i*x;
        V6 = V(6)*exp(1i*delta(6));
        V7 = V(7)*exp(1i*delta(7));
        I  = (V6 - V7) / z;           % 串联电流（pu）
        S_cu = (abs(I)^2 * z) * S_base;     % 铜损 + 漏抗损耗
        % 任务书给定“不计铁损”，故铁损 = 0
        fprintf('\n=== 6→7 变压器损耗 ===\n');
        fprintf('  铜损： P = %.2f MW,  Q = %.2f Mvar\n', real(S_cu), imag(S_cu));
        fprintf('  铁损： P = 0.00 MW,  Q = 0.00 Mvar  （任务书忽略）\n');
        fprintf('  总计： P = %.2f MW,  Q = %.2f Mvar\n', real(S_cu), imag(S_cu));
    end
  
    % 输出总损耗
    fprintf('\n总线路损耗:\t有功: %.2f MW\t无功: %.2f Mvar\n', total_P_loss, total_Q_loss);
end

%% 显示结果函数
function display_results(V, delta, P_calc, Q_calc, P_load, Q_load, S_base)
    fprintf('\n=== 节点电压结果 ===\n');
    fprintf('节点\t电压幅值(pu)\t电压相角(度)\n');
    
    for i = 1:length(V)
        fprintf('%d\t%.4f\t\t%.2f\n', i, V(i), rad2deg(delta(i)));
    end
    
    fprintf('\n=== 节点功率结果 ===\n');
    fprintf('节点\t发电有功(MW)\t发电无功(Mvar)\t负荷有功(MW)\t负荷无功(Mvar)\n');
    
    for i = 1:length(V)
        P_gen = P_calc(i) * S_base + P_load(i) * S_base;
        Q_gen = Q_calc(i) * S_base + Q_load(i) * S_base;
        
        fprintf('%d\t%.2f\t\t%.2f\t\t%.2f\t\t%.2f\n', ...
                i, P_gen, Q_gen, P_load(i)*S_base, Q_load(i)*S_base);
    end
end