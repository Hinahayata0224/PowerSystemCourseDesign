 clc; close all;
%% 系统基准值
Sb = 100; Vb = 220; Zb = Vb^2/Sb;
fprintf('=== 系统基准值 ===\n基准功率: %.0f MVA\n基准电压: %.0f kV\n基准阻抗: %.2f Ω\n\n', Sb, Vb, Zb);
%% 修正的发电机功率计算
P_rated = 3 * 300;
total_load = 230 + 350 + 250 + 400;
P_available = P_rated * 0.90;
P_net = P_available * (1 - 0.07);
power_deficit = total_load - P_net;
fprintf('=== 修正的发电机功率计算 ===\n总负荷需求: %.1f MW\n总装机容量: %.1f MW\n上网功率: %.1f MW\n功率缺额: %.1f MW\n\n', total_load, P_rated, P_net, power_deficit);
%% 负荷水平定义
load_cases = {'正常负荷', '高峰负荷', '低谷负荷'};
load_factors = [1.0, 1.2, 0.7];
results = cell(1, 3);
%% 线路参数计算
r1 = 0.054*0.5; x1 = 0.308*0.5; C1 = 2*0.0116e-6;
f = 50; b1 = 2 * pi * f * C1;
%% 线路数据
linedata = [
    1, 2, 45*r1, 45*x1, 45*b1, 45;
    2, 4,  6*r1,  6*x1,  6*b1,  6;  
    4, 5,  5*r1,  5*x1,  5*b1,  5;
    5, 6, 24*r1, 24*x1, 24*b1, 24;
    6, 3, 18*r1, 18*x1, 18*b1, 18;
    3, 2, 17*r1, 17*x1, 17*b1, 17;
];
nline = size(linedata,1);
%% 转换为标幺值
for i = 1:nline
    linedata(i,3) = linedata(i,3)/Zb;
    linedata(i,4) = linedata(i,4)/Zb;
    linedata(i,5) = linedata(i,5)*Zb;
end
%% 形成节点导纳矩阵Ybus（循环外计算，不变）
Ybus = zeros(6,6);
for k = 1:nline
    from = linedata(k,1); to = linedata(k,2);
    r = linedata(k,3); x = linedata(k,4); b = linedata(k,5);
    z = r + 1j*x; y = 1/z;
    Ybus(from,to) = Ybus(from,to) - y;
    Ybus(to,from) = Ybus(to,from) - y;
    Ybus(from,from) = Ybus(from,from) + y + 1j*b/2;
    Ybus(to,to) = Ybus(to,to) + y + 1j*b/2;
end
G = real(Ybus); B = imag(Ybus);

%% 循环计算不同负荷情况
for case_idx = 1:3
    fprintf('\n=== 正在进行 %s 潮流计算 ===\n', load_cases{case_idx});
    load_factor = load_factors(case_idx); 
    %% 节点数据
    busdata = [
        1, 3, 1.00, 0,   0,   0,   0,   0, 1.05, 0.95;
        2, 1, 1.00, 0, 230*load_factor, 120*load_factor,   0,   0, 1.05, 0.95;
        3, 1, 1.00, 0, 350*load_factor, 190*load_factor,   0,   0, 1.05, 0.95;  
        4, 1, 1.00, 0, 250*load_factor,  85*load_factor,   0,   0, 1.05, 0.95;
        5, 1, 1.00, 0, 400*load_factor, 195*load_factor,   0,   0, 1.05, 0.95;
        6, 2, 1.05, 0,   0,   0, P_net, 0, 1.05, 0.95;
    ];
    nbus = size(busdata,1);
    busdata(:,5:8) = busdata(:,5:8)/Sb;
    %% 牛顿-拉夫逊法潮流计算
    max_iter = 500; tolerance = 1e-6;
    V = busdata(:,3); theta = busdata(:,4);
    pq_bus = find(busdata(:,2) == 1); V(pq_bus) = 0.98;
    pv_bus = find(busdata(:,2) == 2); ref_bus = find(busdata(:,2) == 3);
    npq = length(pq_bus); npv = length(pv_bus);

    iter = 0; converged = false; error_history = [];
    while iter < max_iter && ~converged
        iter = iter + 1;
        P_calc = zeros(nbus,1); Q_calc = zeros(nbus,1);
        for i = 1:nbus
            for j = 1:nbus
                P_calc(i) = P_calc(i) + V(i)*V(j)*(G(i,j)*cos(theta(i)-theta(j)) + B(i,j)*sin(theta(i)-theta(j)));
                Q_calc(i) = Q_calc(i) + V(i)*V(j)*(G(i,j)*sin(theta(i)-theta(j)) - B(i,j)*cos(theta(i)-theta(j)));
            end
        end
        dP = busdata(:,7) - busdata(:,5) - P_calc;
        dQ = busdata(:,8) - busdata(:,6) - Q_calc;
        dQ(pv_bus) = 0; dQ(ref_bus) = 0; dP(ref_bus) = 0;
        max_error = max([abs(dP); abs(dQ)]); error_history(iter) = max_error;
        if max_error < tolerance, converged = true; break; end

        % 形成雅可比矩阵
        J11 = zeros(nbus-1, nbus-1); J12 = zeros(nbus-1, npq);
        J21 = zeros(npq, nbus-1); J22 = zeros(npq, npq);
        
        for i = 1:nbus-1
            m = i + (i >= ref_bus);
            for j = 1:nbus-1
                n = j + (j >= ref_bus);
                if m == n
                    J11(i,j) = -Q_calc(m) - V(m)^2*B(m,m);
                else
                    J11(i,j) = V(m)*V(n)*(G(m,n)*sin(theta(m)-theta(n)) - B(m,n)*cos(theta(m)-theta(n)));
                end
            end
        end
        
        for i = 1:nbus-1
            m = i + (i >= ref_bus);
            for j = 1:npq
                n = pq_bus(j);
                if m == n
                    J12(i,j) = P_calc(m)/V(m) + V(m)*G(m,m);
                else
                    J12(i,j) = V(m)*(G(m,n)*cos(theta(m)-theta(n)) + B(m,n)*sin(theta(m)-theta(n)));
                end
            end
        end
        
        for i = 1:npq
            m = pq_bus(i);
            for j = 1:nbus-1
                n = j + (j >= ref_bus);
                if m == n
                    J21(i,j) = P_calc(m) - V(m)^2*G(m,m);
                else
                    J21(i,j) = -V(m)*V(n)*(G(m,n)*cos(theta(m)-theta(n)) + B(m,n)*sin(theta(m)-theta(n)));
                end
            end
        end
        
        for i = 1:npq
            m = pq_bus(i);
            for j = 1:npq
                n = pq_bus(j);
                if m == n
                    J22(i,j) = Q_calc(m)/V(m) - V(m)*B(m,m);
                else
                    J22(i,j) = V(m)*(G(m,n)*sin(theta(m)-theta(n)) - B(m,n)*cos(theta(m)-theta(n)));
                end
            end
        end

        J = [J11, J12; J21, J22];
        mismatch = [dP(1:ref_bus-1); dP(ref_bus+1:end); dQ(pq_bus)];
        correction = J \ mismatch;

        dtheta = correction(1:nbus-1);
        theta(1:ref_bus-1) = theta(1:ref_bus-1) + dtheta(1:ref_bus-1);
        theta(ref_bus+1:end) = theta(ref_bus+1:end) + dtheta(ref_bus:end);
        dV = correction(nbus:end); V(pq_bus) = V(pq_bus) + dV;
        V(pq_bus) = min(max(V(pq_bus), busdata(pq_bus,10)), busdata(pq_bus,9));
    end

    %% 存储结果
    case_result.converged = converged; case_result.iter = iter;
    case_result.V = V; case_result.theta = theta; case_result.busdata = busdata;
    case_result.load_factor = load_factor;
    
    if converged
        % 计算最终功率
        P_final = zeros(nbus,1); Q_final = zeros(nbus,1);
        for i = 1:nbus
            for j = 1:nbus
                P_final(i) = P_final(i) + V(i)*V(j)*(G(i,j)*cos(theta(i)-theta(j)) + B(i,j)*sin(theta(i)-theta(j)));
                Q_final(i) = Q_final(i) + V(i)*V(j)*(G(i,j)*sin(theta(i)-theta(j)) - B(i,j)*cos(theta(i)-theta(j)));
            end
        end
        
        P_final_MW = P_final * Sb; Q_final_Mvar = Q_final * Sb;
        case_result.P_final = P_final; case_result.Q_final = Q_final;
        case_result.P_final_MW = P_final_MW; case_result.Q_final_Mvar = Q_final_Mvar;
        %% 输出结果
        fprintf('✓ %s潮流计算收敛于第%d次迭代\n', load_cases{case_idx}, iter);
        fprintf('节点\t类型\t电压(pu)\t角度(度)\t负荷(MW)\t负荷(Mvar)\t发电(MW)\t发电(Mvar)\n');
        for i = 1:nbus
            type_str = ''; switch busdata(i,2), case 1, type_str = 'PQ'; case 2, type_str = 'PV'; case 3, type_str = '平衡'; end
            P_load = busdata(i,5) * Sb; Q_load = busdata(i,6) * Sb;
            P_gen = P_final(i) * Sb + P_load; Q_gen = Q_final(i) * Sb + Q_load;
            fprintf('%d\t%s\t%.4f\t\t%6.2f\t\t%6.1f\t\t%6.1f\t\t%6.1f\t\t%6.1f\n', i, type_str, V(i), theta(i)*180/pi, P_load, Q_load, P_gen, Q_gen);
        end

        %% 功率平衡检查
        total_load_P = sum(busdata(:,5)) * Sb; total_load_Q = sum(busdata(:,6)) * Sb;
        total_gen_P = 0; total_gen_Q = 0;
        for i = 1:nbus
            total_gen_P = total_gen_P + (P_final(i) * Sb + busdata(i,5) * Sb);
            total_gen_Q = total_gen_Q + (Q_final(i) * Sb + busdata(i,6) * Sb);
        end
        total_loss_P = total_gen_P - total_load_P; total_loss_Q = total_gen_Q - total_load_Q;
        fprintf('总发电: P=%.2f MW, Q=%.2f Mvar\n总负荷: P=%.2f MW, Q=%.2f Mvar\n网损: P=%.2f MW, Q=%.2f Mvar\n网损率: P=%.2f%%, Q=%.2f%%\n', total_gen_P, total_gen_Q, total_load_P, total_load_Q, total_loss_P, total_loss_Q, total_loss_P/total_gen_P*100, total_loss_Q/total_gen_Q*100);
        
        case_result.total_load_P = total_load_P; case_result.total_load_Q = total_load_Q;
        case_result.total_gen_P = total_gen_P; case_result.total_gen_Q = total_gen_Q;
        case_result.total_loss_P = total_loss_P; case_result.total_loss_Q = total_loss_Q;
        %% 线路潮流分析
        fprintf('线路\t\t有功(MW)\t无功(Mvar)\t电流(A)\t\t负载率(%%)\n');
        total_line_loss_P = 0; line_results = [];
        for k = 1:nline
            from = linedata(k,1); to = linedata(k,2);
            P_flow_from = V(from)^2*G(from,to) - V(from)*V(to)*(G(from,to)*cos(theta(from)-theta(to)) + B(from,to)*sin(theta(from)-theta(to)));
            Q_flow_from = -V(from)^2*B(from,to) - V(from)*V(to)*(G(from,to)*sin(theta(from)-theta(to)) - B(from,to)*cos(theta(from)-theta(to)));
            I_flow = sqrt(P_flow_from^2 + Q_flow_from^2) * Sb / (sqrt(3) * Vb * V(from)) * 1000;
            I_max = 1190; loading_percent = I_flow / I_max * 100;
            fprintf('%d-%d\t\t%7.2f\t\t%7.2f\t\t%7.1f\t\t%6.1f\n', from, to, P_flow_from*Sb, Q_flow_from*Sb, I_flow, loading_percent);
            if loading_percent > 100, fprintf('⚠️ 警告: 线路 %d-%d 过载! (%.1f%%)\n', from, to, loading_percent); end
            P_flow_to = V(to)^2*G(to,from) - V(to)*V(from)*(G(to,from)*cos(theta(to)-theta(from)) + B(to,from)*sin(theta(to)-theta(from)));
            total_line_loss_P = total_line_loss_P + (P_flow_from + P_flow_to) * Sb;
            line_results(k).from = from; line_results(k).to = to; line_results(k).P_flow = P_flow_from*Sb;
            line_results(k).Q_flow = Q_flow_from*Sb; line_results(k).I_flow = I_flow; line_results(k).loading = loading_percent;
        end
        fprintf('线路总有功损耗: %.2f MW\n', total_line_loss_P);
        case_result.line_results = line_results; case_result.total_line_loss_P = total_line_loss_P;
        %% 电压水平分析
        fprintf('节点电压水平分析:\n'); voltage_results = [];
        for i = 1:nbus
            V_kV = V(i) * Vb; status = '正常';
            if V(i) < busdata(i,10), status = '⚠️ 电压偏低'; elseif V(i) > busdata(i,9), status = '⚠️ 电压偏高'; end
            fprintf('节点 %d: %.2f kV (%.3f pu) - %s\n', i, V_kV, V(i), status);
            voltage_results(i).node = i; voltage_results(i).V_pu = V(i); voltage_results(i).V_kV = V_kV; voltage_results(i).status = status;
        end
        case_result.voltage_results = voltage_results;
        %% 系统运行状态评估
        fprintf('发电能力: %.1f MW / %.1f MW (%.1f%%)\n', P_net, P_rated, P_net/P_rated*100);
        fprintf('负荷需求: %.1f MW\n', total_load_P);
        fprintf('功率缺额: %.1f MW\n', total_load_P - P_net);
        fprintf('500kV系统注入功率: %.1f MW\n', P_final(ref_bus)*Sb);
        if total_load_P > P_net, fprintf('🔴 系统状态: 发电不足\n'); else, fprintf('🟢 系统状态: 发电充足\n'); end
        case_result.power_deficit = total_load_P - P_net; case_result.injected_power = P_final(ref_bus)*Sb;
    else
        fprintf('✗ %s潮流计算在%d次迭代后未收敛\n', load_cases{case_idx}, max_iter);
    end
    results{case_idx} = case_result;
end
%% 不同负荷情况对比分析
fprintf('\n=== 不同负荷情况综合对比分析 ===\n');
fprintf('1. 节点电压对比分析:\n节点\t正常负荷(pu)\t高峰负荷(pu)\t低谷负荷(pu)\n');
for i = 2:6
    V_normal = results{1}.V(i); V_peak = results{2}.V(i); V_valley = results{3}.V(i);
    fprintf('%d\t%.4f\t\t%.4f\t\t%.4f\n', i, V_normal, V_peak, V_valley);
end

fprintf('\n2. 功率平衡对比分析:\n负荷情况\t总负荷(MW)\t发电功率(MW)\t功率缺额(MW)\n');
for i = 1:3
    if results{i}.converged
        fprintf('%s\t%.1f\t\t%.1f\t\t%.1f\n', load_cases{i}, results{i}.total_load_P, P_net, results{i}.power_deficit);
    end
end

fprintf('\n3. 系统网损对比分析:\n负荷情况\t总有功网损(MW)\t有功网损率(%%)\n');
for i = 1:3
    if results{i}.converged
        loss_rate = results{i}.total_loss_P / results{i}.total_gen_P * 100;
        fprintf('%s\t%.2f\t\t\t%.2f\n', load_cases{i}, results{i}.total_loss_P, loss_rate);
    end
end

fprintf('\n4. 关键线路负载率分析:\n线路\t\t正常负荷(%%)\t高峰负荷(%%)\t低谷负荷(%%)\n');
critical_lines = [1, 2, 4, 6];
for k = critical_lines
    from = linedata(k,1); to = linedata(k,2);
    loadings = [];
    for i = 1:3
        if results{i}.converged, loadings(i) = results{i}.line_results(k).loading; else, loadings(i) = NaN; end
    end
    fprintf('%d-%d\t\t%.1f\t\t%.1f\t\t%.1f\n', from, to, loadings(1), loadings(2), loadings(3));
end
%% 运行建议
fprintf('\n5. 系统运行建议:\n');
if results{2}.converged
    low_voltage_nodes = []; overloaded_lines = [];
    for i = 1:length(results{2}.voltage_results)
        if results{2}.voltage_results(i).V_pu < 0.95 && i ~= 1, low_voltage_nodes = [low_voltage_nodes, i]; end
    end
    for k = 1:length(results{2}.line_results)
        if results{2}.line_results(k).loading > 100, overloaded_lines = [overloaded_lines, k]; end
    end
    if ~isempty(low_voltage_nodes), fprintf('• 高峰负荷时电压偏低节点: '); fprintf('%d ', low_voltage_nodes); fprintf('\n'); end
    if ~isempty(overloaded_lines), fprintf('• 高峰负荷时过载线路: '); for k = overloaded_lines, fprintf('%d-%d ', linedata(k,1), linedata(k,2)); end; fprintf('\n'); end
    if results{2}.power_deficit > 0.1 * P_net, fprintf('• 高峰负荷时功率缺额较大 (%.1f MW)\n', results{2}.power_deficit); end
end
%% 绘图函数
if results{1}.converged
    plot_power_flow_complete(results{1}.busdata, linedata, results{1}.V, results{1}.theta, Sb, Vb, '正常负荷');
    plot_comparison_analysis(results, load_cases, linedata);
end
%% 潮流图绘制函数
function plot_power_flow_complete(busdata, linedata, V, theta, Sb, Vb, case_name)
    figure('Position', [50, 30, 1300, 500]);
    node_pos = [1, 0.1, 0.5; 2, 0.3, 0.7; 3, 0.3, 0.3; 4, 0.5, 0.7; 5, 0.7, 0.7; 6, 0.7, 0.5;];
    hold on;
    for i = 1:size(linedata,1)
        from = linedata(i,1); to = linedata(i,2);
        x1 = node_pos(from,2); y1 = node_pos(from,3); x2 = node_pos(to,2); y2 = node_pos(to,3);
        P_flow = calculate_line_flow(from, to, V, theta, linedata, i) * Sb;
        Q_flow = calculate_line_flow_q(from, to, V, theta, linedata, i) * Sb;
        line_width = 1 + abs(P_flow)/300;
        color_intensity = min(1, abs(P_flow)/500);
        if P_flow > 0, line_color = [color_intensity, 0, 0]; else, line_color = [0, 0, color_intensity]; end
        plot([x1, x2], [y1, y2], '-', 'LineWidth', line_width, 'Color', line_color);
        mid_x = (x1+x2)/2; mid_y = (y1+y2)/2;
        if P_flow >= 0
            text(mid_x, mid_y, sprintf('→%.1fMW\n%.1fMvar', P_flow, Q_flow), 'FontSize', 8, 'BackgroundColor', 'white', 'HorizontalAlignment', 'center');
        else
            text(mid_x, mid_y, sprintf('←%.1fMW\n%.1fMvar', -P_flow, -Q_flow), 'FontSize', 8, 'BackgroundColor', 'white', 'HorizontalAlignment', 'center');
        end
    end
    for i = 1:size(busdata,1)
        x = node_pos(i,2); y = node_pos(i,3);
        switch busdata(i,2)
            case 1, color = 'blue'; marker = 'o'; node_type = 'PQ';
            case 2, color = 'red'; marker = 's'; node_type = 'PV';  
            case 3, color = 'green'; marker = '^'; node_type = '平衡';
        end
        plot(x, y, marker, 'MarkerSize', 16, 'MarkerFaceColor', color, 'MarkerEdgeColor', 'k', 'LineWidth', 2);
        V_kV = V(i) * Vb;
        text(x, y+0.08, sprintf('%s节点 %d\nV=%.3f pu\n%.1f kV∠%.1f°', node_type, i, V(i), V_kV, theta(i)*180/pi), 'HorizontalAlignment', 'center', 'FontSize', 9, 'BackgroundColor', 'white', 'EdgeColor', 'black');
        P_load = busdata(i,5) * Sb; Q_load = busdata(i,6) * Sb; P_gen = busdata(i,7) * Sb;
        if P_load > 0 || Q_load > 0
            text(x, y-0.08, sprintf('负荷: %.0f+j%.0f', P_load, Q_load), 'HorizontalAlignment', 'center', 'FontSize', 8, 'BackgroundColor', 'white', 'EdgeColor', 'red');
        end
        if P_gen > 0
            text(x, y-0.12, sprintf('发电: %.0fMW', P_gen), 'HorizontalAlignment', 'center', 'FontSize', 8, 'BackgroundColor', 'white', 'EdgeColor', 'green');
        end
    end
    title(['电力系统潮流分布图 - ' case_name '情况'], 'FontSize', 14, 'FontWeight', 'bold');
    legend_elements = [plot(NaN, NaN, 'o', 'MarkerSize', 10, 'MarkerFaceColor', 'blue', 'MarkerEdgeColor', 'k'), plot(NaN, NaN, 's', 'MarkerSize', 10, 'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'k'), plot(NaN, NaN, '^', 'MarkerSize', 10, 'MarkerFaceColor', 'green', 'MarkerEdgeColor', 'k')];
    legend(legend_elements, {'PQ节点', 'PV节点', '平衡节点'}, 'Location', 'northeast');
    axis equal; xlim([0, 0.9]); ylim([0.2, 0.8]); grid on;
    text(0.02, 0.15, sprintf('系统基准值:\nS_b = %.0f MVA\nV_b = %.0f kV\n\n线路功率: MW + jMvar\n箭头表示功率流向', Sb, Vb), 'FontSize', 10, 'BackgroundColor', 'white', 'EdgeColor', 'black', 'VerticalAlignment', 'bottom');
end

function P_flow = calculate_line_flow(from, to, V, theta, linedata, line_idx)
    r = linedata(line_idx,3); x = linedata(line_idx,4); z = r + 1j*x; y = 1/z;
    P_flow = V(from)^2*real(y) - V(from)*V(to)*(real(y)*cos(theta(from)-theta(to)) + imag(y)*sin(theta(from)-theta(to)));
end

function Q_flow = calculate_line_flow_q(from, to, V, theta, linedata, line_idx)
    r = linedata(line_idx,3); x = linedata(line_idx,4); b = linedata(line_idx,5); z = r + 1j*x; y = 1/z;
    Q_flow = -V(from)^2*(imag(y) + b/2) - V(from)*V(to)*(real(y)*sin(theta(from)-theta(to)) - imag(y)*cos(theta(from)-theta(to)));
end
%% 综合对比分析绘图函数
function plot_comparison_analysis(results, load_cases, linedata)
    figure('Position', [100, 100, 1400, 900]);
    
    subplot(2, 3, 1);
    nodes = 1:6; voltages = zeros(3, 6);
    for i = 1:3, if results{i}.converged, voltages(i, :) = results{i}.V'; end, end
    plot(nodes, voltages(1, :), 'bo-', 'LineWidth', 2, 'MarkerSize', 8); hold on;
    plot(nodes, voltages(2, :), 'rs-', 'LineWidth', 2, 'MarkerSize', 8);
    plot(nodes, voltages(3, :), 'g^-', 'LineWidth', 2, 'MarkerSize', 8);
    grid on; xlabel('节点编号'); ylabel('电压 (pu)'); title('节点电压对比');
    legend(load_cases, 'Location', 'best'); set(gca, 'XTick', nodes);
    ylim([0.94, 1.06]); plot([0.5, 6.5], [0.95, 0.95], 'r--', 'LineWidth', 1); plot([0.5, 6.5], [1.05, 1.05], 'r--', 'LineWidth', 1);
    
    subplot(2, 3, 2);
    power_data = zeros(1, 3); for i = 1:3, if results{i}.converged, power_data(i) = results{i}.power_deficit; end, end
    bar(power_data, 'FaceColor', [0.7, 0.7, 0.9]); set(gca, 'XTickLabel', load_cases); ylabel('功率缺额 (MW)'); title('功率缺额对比'); grid on;
    for i = 1:3, text(i, power_data(i) + max(power_data)*0.05, sprintf('%.1f', power_data(i)), 'HorizontalAlignment', 'center', 'FontWeight', 'bold'); end
    
    subplot(2, 3, 3);
    loss_data = zeros(1, 3); for i = 1:3, if results{i}.converged, loss_data(i) = results{i}.total_loss_P; end, end
    bar(loss_data, 'FaceColor', [0.9, 0.7, 0.7]); set(gca, 'XTickLabel', load_cases); ylabel('总有功网损 (MW)'); title('网损对比'); grid on;
    for i = 1:3, text(i, loss_data(i) + max(loss_data)*0.05, sprintf('%.1f', loss_data(i)), 'HorizontalAlignment', 'center', 'FontWeight', 'bold'); end
    
    subplot(2, 3, 4);
    critical_lines = [1, 2, 4, 6]; line_labels = {}; loading_data = zeros(3, length(critical_lines));
    for j = 1:length(critical_lines)
        k = critical_lines(j); line_labels{j} = sprintf('%d-%d', linedata(k,1), linedata(k,2));
        for i = 1:3, if results{i}.converged, loading_data(i, j) = results{i}.line_results(k).loading; end, end
    end
    x = 1:length(critical_lines); bar_handles = bar(x, loading_data');
    set(gca, 'XTickLabel', line_labels); ylabel('负载率 (%)'); title('关键线路负载率对比'); legend(load_cases, 'Location', 'northwest'); grid on;
    ylim([0, max(max(loading_data))*1.1]); plot([0, length(critical_lines)+1], [100, 100], 'r--', 'LineWidth', 2);
    
    subplot(2, 3, 5);
    injected_power = zeros(1, 3); for i = 1:3, if results{i}.converged, injected_power(i) = results{i}.injected_power; end, end
    bar(injected_power, 'FaceColor', [0.7, 0.9, 0.7]); set(gca, 'XTickLabel', load_cases); ylabel('500kV注入功率 (MW)'); title('500kV注入功率对比'); grid on;
    for i = 1:3, text(i, injected_power(i) + max(injected_power)*0.05, sprintf('%.1f', injected_power(i)), 'HorizontalAlignment', 'center', 'FontWeight', 'bold'); end
    
    subplot(2, 3, 6); axis off;
    summary_text = sprintf('系统运行状态总结:\n\n');
    for i = 1:3
        if results{i}.converged
            low_voltage_count = 0; overload_count = 0;
            for j = 1:length(results{i}.voltage_results), if results{i}.voltage_results(j).V_pu < 0.95 && j ~= 1, low_voltage_count = low_voltage_count + 1; end, end
            for k = 1:length(results{i}.line_results), if results{i}.line_results(k).loading > 100, overload_count = overload_count + 1; end, end
            status = '良好'; if low_voltage_count > 0 || overload_count > 0 || results{i}.power_deficit > 0.1*753.39, status = '需关注'; end
            if low_voltage_count > 1 || overload_count > 1 || results{i}.power_deficit > 0.2*753.39, status = '较差'; end
            summary_text = [summary_text, sprintf('%s:\n• 电压异常节点: %d个\n• 过载线路: %d条\n• 功率缺额: %.1f MW\n• 状态评估: %s\n\n', load_cases{i}, low_voltage_count, overload_count, results{i}.power_deficit, status)];
        else
            summary_text = [summary_text, sprintf('%s: 计算未收敛\n\n', load_cases{i})];
        end
    end
    text(0.1, 0.9, summary_text, 'FontSize', 11, 'VerticalAlignment', 'top', 'BackgroundColor', [0.95, 0.95, 0.95], 'EdgeColor', 'black', 'Margin', 10);
    sgtitle('电力系统不同负荷情况下运行状态综合对比分析', 'FontSize', 16, 'FontWeight', 'bold');
end