function power_system_stability_analysis()
% 电力系统稳定计算程序 - 发电机二阶经典模型

clc;
clear;
close all;

%% 1. 运行潮流计算获取初始状态
fprintf('=== 电力系统稳定分析 ===\n');
fprintf('步骤1: 运行潮流计算获取初始状态...\n');

[V, delta, P_gen, Q_gen, iter] = C2task2_powerflow();
fprintf('潮流计算完成，迭代次数: %d\n', iter);

%% 2. 系统参数设置
fprintf('\n步骤2: 设置系统参数...\n');

S_base = 1230; % MVA
f0 = 50; % Hz
omega0 = 2 * pi * f0;

% 发电机参数（三台机组合并等效）
P_rated = 3 * 300; % 额定功率 900 MW
S_rated = 3 * 350; % 额定容量 1050 MVA
V_rated = 10.5; % kV
PF = 0.85; % 功率因数

% 发电机二阶模型参数（标幺值，基准1230MVA）
Xd = 1.8 * (S_base / 1050); % d轴同步电抗
Xd_prime = 0.18 * (S_base / 1050); % d轴暂态电抗
Xq = 1.2 * (S_base / 1050); % q轴同步电抗
Tj = 7; % 惯性时间常数（秒）
H = Tj / 2; % 惯性常数
D = 1; % 阻尼系数
    
%% 3. 计算发电机初始状态
fprintf('\n步骤3: 计算发电机初始状态...\n');

% 发电机节点（节点7）
gen_bus = 7;
slack_bus = 1; % 平衡节点（无穷大系统）

% 获取发电机端电压和功率
V_gen = V(gen_bus);
delta_gen = delta(gen_bus);
P_gen_actual = P_gen(gen_bus) ; % MW
Q_gen_actual = Q_gen(gen_bus); % Mvar

fprintf('发电机初始状态:\n');
fprintf('  端电压: %.4f pu (%.2f kV)\n', V_gen, V_gen * 10.5);
fprintf('  功角: %.2f 度\n', rad2deg(delta_gen));
fprintf('  出力: P = %.2f MW, Q = %.2f Mvar\n', P_gen_actual, Q_gen_actual);

% 计算发电机内电势（经典模型）
[E_prime, delta0] = calculate_internal_voltage(V_gen, delta_gen, P_gen_actual/S_base, ...
    Q_gen_actual/S_base, Xd_prime);

fprintf('发电机内电势（经典模型）:\n');
fprintf('  E'' = %.4f pu, δ0 = %.2f 度\n', E_prime, rad2deg(delta0));


%% 4. 网络变换法计算转移阻抗
fprintf('\n步骤4: 网络变换法计算转移阻抗...\n');

% 正常情况转移阻抗
X_transfer = calculate_transfer_impedance_star_delta();
fprintf('  故障前转移阻抗:X_transfer = %.4f pu\n', X_transfer);

% 计算故障期间的转移阻抗（节点2三相短路）
X_transfer_fault = calculate_transfer_impedance_fault_star_delta();
fprintf('  故障期间转移阻抗:X_transfer_fault = %.6f pu\n', X_transfer_fault);

% 计算故障后转移阻抗（切除线路3-6的其中一回）
X_transfer_post = calculate_transfer_impedance_post_star_delta();
fprintf('  故障后转移阻抗:X_transfer_post = %.6f pu\n', X_transfer_post);

%% 5. 稳定计算
fprintf('\n步骤5: 进行稳定计算...\n');

% 故障设置（三相短路）- 修改为节点2
fault_bus = 2; % 故障发生在节点2（修改点）
fault_start = 0.5; % 故障开始时间
fault_duration = 0.15; % 故障持续时间
fault_cleared = fault_start + fault_duration; % 故障切除时间

fprintf('  故障设置:三相短路\n');
fprintf('  故障节点: %d\n', fault_bus);
fprintf('  故障开始: %.2f s, 切除: %.2f s\n', fault_start, fault_cleared);

% 仿真参数
t_start = 0;
t_end = 20; % 仿真时间 20秒
t_step = 0.01; % 时间步长

% 初始状态向量 [功角(rad), 角速度偏差(pu)]
x0 = [delta0; 0]; % 初始角速度偏差为0
P_m0 = P_gen_actual / S_base;

% 进行稳定计算
[t, x] = run_stability_simulation(x0, t_start, t_end, t_step, E_prime, ...
    X_transfer, X_transfer_fault, X_transfer_post, H, D, omega0, fault_start, fault_cleared, P_m0);

%% 6. 结果输出和绘图
fprintf('\n步骤6: 输出结果和绘制曲线...\n');

% 提取结果
delta_t = x(:, 1); % 功角轨迹 (rad)
omega_t = x(:, 2); % 角速度偏差 (pu)
omega_actual = omega_t + 1; % 实际角速度 (pu)

% 转换为实际值
delta_deg = rad2deg(delta_t); % 功角 (度)
omega_rpm = omega_actual * 3000; % 转速 (rpm) - 对于50Hz系统

% 显示关键结果
fprintf('\n=== 稳定计算结果 ===\n');
fprintf('最大功角: %.2f 度 (在 t=%.2f s)\n', max(delta_deg), t(delta_deg == max(delta_deg)));
fprintf('最小功角: %.2f 度 (在 t=%.2f s)\n', min(delta_deg), t(delta_deg == min(delta_deg)));
fprintf('最大转速: %.2f rpm\n', max(omega_rpm));
fprintf('最小转速: %.2f rpm\n', min(omega_rpm));

% 绘制功角和角速度曲线
plot_stability_results(t, delta_deg, omega_rpm, fault_start, fault_cleared);

%% 7. 计算极限切除时间和极限切除角
fprintf('\n步骤7: 计算极限切除时间和极限切除角...\n');

% 计算机械功率
P_m = P_m0; % 初始机械功率
fprintf('初始机械功率: %.4f pu\n', P_m);

% 通过逐步增加故障切除时间来确定极限切除时间
fault_duration_trial = 0.01:0.005:0.3; % 故障持续时间试验值
n_trials = length(fault_duration_trial);
stability_flags = false(n_trials, 1); % 稳定性标志

fprintf('进行极限切除时间搜索...\n');
for i = 1:n_trials
    fault_duration_current = fault_duration_trial(i);
    fault_cleared_current = fault_start + fault_duration_current;
    
    % 进行稳定计算
    [t_current, x_current] = run_stability_simulation(x0, t_start, t_end, t_step, E_prime, ...
        X_transfer, X_transfer_fault, X_transfer_post,H, D, omega0, fault_start, fault_cleared_current,P_m0);
    
    % 判断稳定性（功角是否收敛）
    stable = is_stable(x_current(:,1), t_current);
    stability_flags(i) = stable;
 
end
fprintf('极限切除时间搜索完成\n');

% 找到极限切除时间
critical_indices = find(stability_flags);
if ~isempty(critical_indices)
    critical_clearing_time = fault_duration_trial(critical_indices(end));
    critical_clearing_angle = calculate_critical_clearing_angle(E_prime, X_transfer, X_transfer_fault, X_transfer_post, P_m, delta0);
    
    fprintf('\n=== 极限切除分析结果 ===\n');
    fprintf('极限切除时间: %.3f s\n', critical_clearing_time);
    fprintf('极限切除角: %.2f 度\n', rad2deg(critical_clearing_angle));
else
    fprintf('警告：未找到稳定的故障切除时间！\n');
end

fprintf('\n稳定计算完成！\n');

end

%% 辅助函数定义
% 计算发电机内电势（经典二阶模型）
function [E_prime, delta_int] = calculate_internal_voltage(Vt, delta_t, P, Q, Xd_prime)
% 输入: Vt - 端电压, delta_t - 功角, P - 有功, Q - 无功, Xd_prime - d轴暂态电抗
% 输出: E_prime - 内电势幅值, delta_int - 内电势角度

% 计算电流
 I = (P - 1j*Q) / conj(Vt);  % 电流计算修正
    E_prime_complex = Vt + 1j * Xd_prime * I;
    E_prime = abs(E_prime_complex);
    delta_int = angle(E_prime_complex);
end

%% 星三角变换方法计算转移阻抗
function X_transfer = calculate_transfer_impedance_star_delta()
    fprintf('  使用网络变换法计算正常情况转移阻抗...\n');

     % 系统基准值
    S_base = 1230;
    V_base_230 = 230;
    Z_base_230 = (V_base_230^2) / S_base;

    % 线路参数
    r_ohm_km = 0.054;
    x_ohm_km = 0.308;
    r_pu_km = r_ohm_km / Z_base_230;
    x_pu_km = x_ohm_km / Z_base_230;

    % 发电机参数
    Xd_prime = 0.18 * (S_base / 1050);
    X_T = 0.105 * (S_base / 1050);

    % 线路阻抗（忽略电阻，只考虑电抗）
    Z_12 = 1j * x_pu_km * 30;  % 线路1-2
    Z_23 = 1j * x_pu_km * 17;  % 线路2-3
    Z_24 = 1j * x_pu_km * 6;   % 线路2-4
    Z_36 = 1j * x_pu_km * 18;  % 线路3-6
    Z_45 = 1j * x_pu_km * 5;   % 线路4-5
    Z_65 = 1j * x_pu_km * 24;  % 线路6-5

    % 变压器和发电机阻抗
    Z_T = 1j * X_T;            % 变压器6-7
    Z_G = 1j * Xd_prime;       % 发电机7-8

    fprintf('    各支路电抗值 (pu):\n');
    fprintf('      Z12 = j%.4f, Z23 = j%.4f, Z24 = j%.4f\n', imag(Z_12), imag(Z_23), imag(Z_24));
    fprintf('      Z36 = j%.4f, Z45 = j%.4f, Z65 = j%.4f\n', imag(Z_36), imag(Z_45), imag(Z_65));
    fprintf('      ZT = j%.4f, ZG = j%.4f\n', imag(Z_T), imag(Z_G));

    % 网络化简
    Z_path1_26 = Z_23 + Z_36;
    Z_path2_26 = Z_24 + Z_45 + Z_65;
    Z_26_combined = 1 / (1/Z_path1_26 + 1/Z_path2_26);
    fprintf('    节点2到节点6的并联阻抗:\n');
    fprintf('      路径1 (2-3-6): Z = j%.4f\n', imag(Z_path1_26));
    fprintf('      路径2 (2-4-5-6): Z = j%.4f\n', imag(Z_path2_26));
    fprintf('      并联等效: Z_26 = j%.4f\n', imag(Z_26_combined));

    Z_16_via_2 = Z_12 + Z_26_combined;
    Z_18_total = Z_16_via_2 + Z_T + Z_G;
    fprintf('      并联等效: Z_18 = j%.4f\n', imag(Z_18_total));
    X_transfer = abs(imag(Z_18_total));
end

function X_transfer_fault = calculate_transfer_impedance_fault_star_delta()
    fprintf('  使用网络变换法计算故障期间转移阻抗...\n');

    S_base = 1230;
    V_base_230 = 230;
    Z_base_230 = (V_base_230^2) / S_base;

    r_ohm_km = 0.054;
    x_ohm_km = 0.308;
    r_pu_km = r_ohm_km / Z_base_230;
    x_pu_km = x_ohm_km / Z_base_230;

    Xd_prime = 0.18 * (S_base / 1050);
    X_T = 0.105 * (S_base / 1050);

    % 线路阻抗
    Z_12 = 1j * x_pu_km * 30;  % 线路1-2
    Z_23 = 1j * x_pu_km * 17;  % 线路2-3
    Z_24 = 1j * x_pu_km * 6;   % 线路2-4
    Z_36 = 1j * x_pu_km * 18;  % 线路3-6
    Z_45 = 1j * x_pu_km * 5;   % 线路4-5
    Z_65 = 1j * x_pu_km * 24;  % 线路6-5

    % 变压器和发电机阻抗
    Z_T = 1j * X_T;            % 变压器6-7
    Z_G = 1j * Xd_prime;       % 发电机7-8

    fprintf('    各支路电抗值 (pu):\n');
    fprintf('      Z12 = j%.4f, Z23 = j%.4f, Z24 = j%.4f\n', imag(Z_12), imag(Z_23), imag(Z_24));
    fprintf('      Z36 = j%.4f, Z45 = j%.4f, Z65 = j%.4f\n', imag(Z_36), imag(Z_45), imag(Z_65));
    fprintf('      ZT = j%.4f, ZG = j%.4f\n', imag(Z_T), imag(Z_G));

    % 网络化简

    % 步骤1: 从发电机内节点8开始，到变压器节点7
    Z_87 = Z_G;  % 发电机阻抗

    % 步骤2: 节点7到节点6（变压器阻抗）
    Z_76 = Z_T;

    % 步骤3: 节点6到节点2的等效阻抗
    % 路径1: 6-3-2
    Z_63 = Z_36;  % 线路3-6
    Z_32 = Z_23;  % 线路2-3
    Z_path1 = Z_63 + Z_32;

    % 路径2: 6-5-4-2
    Z_65_path = Z_65;  % 线路6-5
    Z_54 = Z_45;       % 线路5-4
    Z_42 = Z_24;       % 线路4-2
    Z_path2 = Z_65_path + Z_54 + Z_42;

    % 路径1和路径2并联
    Z_62_parallel = 1 / (1/Z_path1 + 1/Z_path2);

    fprintf('    节点6到节点2的并联阻抗:\n');
    fprintf('      路径1 (6-3-2): Z = j%.4f\n', imag(Z_path1));
    fprintf('      路径2 (6-5-4-2): Z = j%.4f\n', imag(Z_path2));
    fprintf('      并联等效: Z_62 = j%.4f\n', imag(Z_62_parallel));

    % 步骤4: 总阻抗 = Z87 + Z76 + Z62
    Z_82_total = Z_87 + Z_76 + Z_62_parallel;
    X_transfer_fault = abs(imag(Z_82_total))+1e6;
end

function X_transfer_post = calculate_transfer_impedance_post_star_delta()
    fprintf('  使用网络变换法计算故障后转移阻抗...\n');

    S_base = 1230;
    V_base_230 = 230;
    Z_base_230 = (V_base_230^2) / S_base;

    r_ohm_km = 0.054;
    x_ohm_km = 0.308;
    r_pu_km = r_ohm_km / Z_base_230;
    x_pu_km = x_ohm_km / Z_base_230;

    Xd_prime = 0.18 * (S_base / 1050);
    X_T = 0.105 * (S_base / 1050);

     % 线路阻抗 - 线路2-3改为单回线，阻抗加倍，其他线路保持双回线
    Z_12 = 1j * x_pu_km * 30;  % 线路1-2（单回线）
    Z_23 = 1j * x_pu_km * 17 * 2;  % 线路2-3（单回线，阻抗加倍）
    Z_24 = 1j * x_pu_km * 6;   % 线路2-4（双回线）
    Z_36 = 1j * x_pu_km * 18;  % 线路3-6（双回线）
    Z_45 = 1j * x_pu_km * 5;   % 线路4-5（双回线）
    Z_65 = 1j * x_pu_km * 24;  % 线路6-5（双回线）

    % 变压器和发电机阻抗
    Z_T = 1j * X_T;            % 变压器6-7
    Z_G = 1j * Xd_prime;       % 发电机7-8

   % 网络化简
    % 步骤1: 计算节点2到节点6的等效阻抗
    Z_path1_26 = Z_23 + Z_36;
    % 路径2: 2-4-5-6
    Z_path2_26 = Z_24 + Z_45 + Z_65;
    % 两条路径并联
    Z_26_combined = 1 / (1/Z_path1_26 + 1/Z_path2_26);

    fprintf('    节点2到节点6的并联阻抗:\n');
    fprintf('      路径1 (2-3-6): Z = j%.4f\n', imag(Z_path1_26));
    fprintf('      路径2 (2-4-5-6): Z = j%.4f\n', imag(Z_path2_26));
    fprintf('      并联等效: Z_26 = j%.4f\n', imag(Z_26_combined));

    % 步骤2: 现在网络简化为: 1-2, 2-6(等效), 6-7, 7-8
    % 计算节点1到节点6的阻抗
    Z_16_via_2 = Z_12 + Z_26_combined;

    % 步骤3: 节点1到节点8的总阻抗
    Z_18_total = Z_16_via_2 + Z_T + Z_G;
    X_transfer_post = abs(imag(Z_18_total));
end

function [t, x] = run_stability_simulation(x0, t_start, t_end, t_step, E_prime, ...
    X_transfer,X_transfer_fault, X_transfer_post, H, D, omega0, fault_start, fault_cleared,P_m0)
% 进行稳定仿真计算

% 时间向量
t = t_start:t_step:t_end;
n_steps = length(t);

% 初始化状态变量
x = zeros(n_steps, 2);
x(1, :) = x0';

% 使用改进欧拉法进行数值积分
for k = 1:n_steps-1
    current_time = t(k);
    
    % 判断当前网络状态
    if current_time < fault_start
        % 故障前状态
        X_eff = X_transfer;
    elseif current_time < fault_cleared
        % 故障期间
        X_eff =X_transfer_fault; 
    else
        % 故障后状态
        X_eff =X_transfer_post ;
    end
    
    % 计算当前电磁功率
    delta_k = x(k, 1);
    P_e = (E_prime * 1.0 / X_eff) * sin(delta_k); % 无穷大系统电压为1.0pu
    
    % 机械功率（恒定）
    P_m = P_m0;% 直接从潮流结果获取机械功率
    
    % 状态方程：转子运动方程
    omega_k = x(k, 2) + 1; % 实际角速度(pu)
    
    % 状态导数
    dx1 = omega0 * (omega_k - 1);
    dx2 = (P_m - P_e - D * (omega_k - 1)) / (2 * H);
    
    % 改进欧拉法预测-校正
    % 预测步
    x_pred = x(k, :) + [dx1, dx2] * t_step;
    
    % 校正步（使用预测值重新计算导数）
    delta_pred = x_pred(1);
    P_e_pred = (E_prime * 1.0 / X_eff) * sin(delta_pred);
    omega_pred = x_pred(2) + 1;
    
    dx1_pred = omega0 * (omega_pred - 1);
    dx2_pred = (P_m - P_e_pred - D * (omega_pred - 1)) / (2 * H);
    
    % 平均值
    dx_avg = 0.5 * ([dx1, dx2] + [dx1_pred, dx2_pred]);
    
    % 更新状态
    x(k+1, :) = x(k, :) + dx_avg * t_step;
end
% fprintf('  时域仿真完成，共计算 %d 个时间点\n', n_steps);
end

% 绘制功角曲线和转速曲线
function plot_stability_results(t, delta_deg, omega_rpm, fault_start, fault_cleared)
    figure('Position', [100, 100, 1200, 800]);
    
    % 功角曲线
    subplot(2, 1, 1);
    plot(t, delta_deg, 'b-', 'LineWidth', 2);
    hold on;
    x_fill = [fault_start, fault_cleared, fault_cleared, fault_start];
    y_fill = [min(delta_deg)-10, min(delta_deg)-10, max(delta_deg)+10, max(delta_deg)+10];
    fill(x_fill, y_fill, 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    grid on;
    xlabel('时间 (s)');
    ylabel('功角 (度)');
    title('发电机功角摇摆曲线');
    legend('功角', '故障期间', 'Location', 'best');
    xlim([0, max(t)]);
    
    % 转速曲线
    subplot(2, 1, 2);
    plot(t, omega_rpm, 'r-', 'LineWidth', 2);
    hold on;
    plot([0, max(t)], [3000, 3000], 'k--', 'LineWidth', 1); % 同步转速参考线
    % 标记故障期间
    x_fill = [fault_start, fault_cleared, fault_cleared, fault_start];
    y_fill = [min(omega_rpm)-10, min(omega_rpm)-10, max(omega_rpm)+10, max(omega_rpm)+10];
    fill(x_fill, y_fill, 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    grid on;
    xlabel('时间 (s)');
    ylabel('转速 (rpm)');
    title('发电机转速变化曲线');
    legend('实际转速', '同步转速(3000rpm)', '故障期间', 'Location', 'best');
    xlim([0, max(t)]);
end

function delta_critical = calculate_critical_clearing_angle(E_prime,X_transfer,X_fault, X_post, P_m, delta0)
% 使用等面积法则计算临界切除角

% 故障前功率最大值
P_max_pre = E_prime * 1.0 / X_transfer;

% 故障期间功率最大值
P_max_fault = E_prime * 1.0 /X_fault;

% 故障后功率最大值
P_max_post = E_prime * 1.0 /X_post;

% 计算临界切除角（等面积法则）
if P_m > P_max_post
    error('机械功率大于故障后最大电磁功率，系统无法稳定！');
end
delta_max = pi - asin(P_m / P_max_post); % 最大可能功角
fprintf('δ_max = %.2f°\n', rad2deg(delta_max));

% 等面积法则计算临界切除角
    % 加速面积 = 减速面积
    % ∫_{δ0}^{δ_c} (P_m - P_e_fault) dδ = ∫_{δ_c}^{δ_max} (P_e_post - P_m) dδ

    % 其中：
    % P_e_fault = P_max_fault * sin(δ)
    % P_e_post = P_max_post * sin(δ)

    % 积分后得到：
    % P_m(δ_c - δ_0) + P_max_fault(cosδ_c - cosδ_0) = P_max_post(cosδ_c - cosδ_max) - P_m(δ_max - δ_c)

    % 整理得：
    % P_m(δ_max - δ_0) = P_max_post(cosδ_c - cosδ_max) - P_max_fault(cosδ_c - cosδ_0)

    % 最终公式：
    % cosδ_c = [P_m(δ_max - δ_0) + P_max_post cosδ_max - P_max_fault cosδ_0] / (P_max_post - P_max_fault)

   if abs(P_max_post - P_max_fault) < 1e-6
        % 如果分母接近0，使用简化公式（当故障期间功率接近0时）
        cos_delta_c = cos(delta_max) + (P_m * (delta_max - delta0)) / P_max_post;
        fprintf('使用简化公式计算临界切除角\n');
   else
        % 标准等面积法公式
        numerator = P_m * (delta_max - delta0) + P_max_post * cos(delta_max) - P_max_fault * cos(delta0);
        denominator = P_max_post - P_max_fault;
        cos_delta_c = numerator / denominator;
        fprintf('使用标准等面积法公式计算临界切除角\n');
   end

    % 确保cos值在有效范围内
   if cos_delta_c > 1
       cos_delta_c = 1;
   elseif cos_delta_c < -1
       cos_delta_c = -1;
   end

delta_critical = acos(cos_delta_c);
fprintf('计算得到的临界切除角: %.2f°\n', rad2deg(delta_critical));

% 验证结果合理性
if delta_critical < delta0
    delta_critical = delta0+ 0.1; % 最小比初始角大一些
elseif delta_critical > delta_max
    delta_critical = delta_max - 0.05; % 最大比极限角小一些
end
end

function stable = is_stable(delta_trajectory, t)
    n = length(delta_trajectory);
    if n < 100
        stable = false;
        return;
    end

    % 取最后5秒的数据进行判断
    t_end = t(end);
    if t_end < 5
        % 如果总仿真时间小于5秒，使用最后1/3数据
        start_idx = round(2*n/3);
    else
        % 取最后5秒数据
        start_idx = find(t >= t_end - 5, 1);
        if isempty(start_idx)
            start_idx = round(0.7 * n);
        end
    end

    final_segment = delta_trajectory(start_idx:end);

    % 计算振荡特性
    peak_to_peak = max(final_segment) - min(final_segment);
    delta_mean = mean(final_segment);

    % 更合理的稳定性条件：
    % 1. 最后阶段振荡幅度 <0.25 rad (约15度)
    % 2. 没有持续的单向漂移
    % 3. 振荡中心相对稳定

    condition1 = peak_to_peak < 0.25;  % 15度

    % 检查是否发散：最后一段与中间段比较
    mid_start = round(0.4 * n);
    mid_end = round(0.6 * n);
    if mid_end > mid_start
        mid_segment = delta_trajectory(mid_start:mid_end);
        mid_amplitude = max(mid_segment) - min(mid_segment);
        condition2 = peak_to_peak <= mid_amplitude * 1.5;  % 允许适度增长
    else
        condition2 = true;
    end

    % 检查是否失步：功角变化是否超过180度
    total_change = abs(delta_trajectory(end) - delta_trajectory(1));
    condition3 = total_change < pi;  % 小于180度

    stable = condition1 && condition2 && condition3;
end