function generate_smooth_load_curve()
    % 时间向量 (24小时，每小时一个点)
    time = 0:1:24;
    
    % 负荷曲线数据点，使用多项式拟合生成平滑曲线
    % 这些点是根据图像手动调整的，以匹配图中的趋势
    load_data = [0.78, 0.72, 0.70, 0.68, 0.67, 0.68, 0.70, 0.75, 0.78, 0.80, ...
                  0.82, 0.81, 0.80, 0.78, 0.76, 0.75, 0.74, 0.76, 0.85, 0.90, ...
                  0.95, 0.93, 0.90, 0.85, 0.78];
    
    % 绘制负荷曲线
    figure('Position', [100, 100, 900, 500]);
    plot(time, load_data, 'b-', 'LineWidth', 2);
    hold on;
    
    % 标记19:00的峰值
    plot(19, 0.9, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    text(19, 0.92, '19:00  0.9 pu', 'FontSize', 10, 'VerticalAlignment', 'bottom');
    
    % 图形美化
    grid on;
    xlabel('时间 (小时)', 'FontSize', 12);
    ylabel('负荷(标幺值)', 'FontSize', 12);
    title('典型日负荷曲线图(夏季)', 'FontSize', 14);
    xlim([0, 24]);
    ylim([0.4, 1.0]);
    set(gca, 'XTick', 0:2:24);
    set(gca, 'FontSize', 11);
    set(gca, 'GridAlpha', 0.2);
    
    % 添加背景色区分时段
    fill([0, 6, 6, 0], [0.4, 0.4, 0.8, 0.8], [0.9, 0.95, 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    fill([6, 18, 18, 6], [0.4, 0.4, 0.8, 0.8], [0.95, 0.98, 0.95], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    fill([18, 24, 24, 18], [0.4, 0.4, 0.8, 0.8], [1, 0.95, 0.9], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    
    text(3, 0.85, '夜间', 'HorizontalAlignment', 'center', 'FontSize', 10);
    text(12, 0.85, '日间', 'HorizontalAlignment', 'center', 'FontSize', 10);
    text(21, 0.85, '晚间', 'HorizontalAlignment', 'center', 'FontSize', 10);
    
    % % % 保存图像
    % saveas(gcf, 'smooth_load_curve.png');
end

% 运行函数
generate_smooth_load_curve();