function fig = plot_scenario(D, savepath)
%PLOT_SCENARIO  Top-down tunnel scenario: BS + true trajectory.
%
%   fig = plot_scenario(D)            draws the figure
%   fig = plot_scenario(D, savepath)  also saves a PNG to savepath
%
% Reproduces the slide-7 style scenario (x = lateral, y = along-tunnel),
% colours the trajectory by active array and marks the array switch.

    x = D.ue_pos(1,:);  y = D.ue_pos(2,:);
    a1 = D.ant_idx == 1; a2 = D.ant_idx == 2;

    fig = figure('Color','w','Position',[100 100 900 520]);

    % approximate tunnel walls from the lateral extent of the path
    halfW = max(abs(x)) + 0.5;
    yl = [min(y)-3, max(y)+3];

    hold on; box on; grid on;
    % tunnel walls (excluded from legend)
    plot([-halfW -halfW], yl, 'k-', 'LineWidth', 1.5, 'HandleVisibility','off');
    plot([ halfW  halfW], yl, 'k-', 'LineWidth', 1.5, 'HandleVisibility','off');

    % trajectory, split by active array
    plot(x(a1), y(a1), '-o', 'Color', [0 0.45 0.74], 'MarkerSize', 3, ...
        'MarkerFaceColor', [0 0.45 0.74], 'DisplayName', 'Array 1 (t=1..59)');
    plot(x(a2), y(a2), '-o', 'Color', [0.85 0.33 0.10], 'MarkerSize', 3, ...
        'MarkerFaceColor', [0.85 0.33 0.10], 'DisplayName', 'Array 2 (t=60..113)');

    % start / switch / end markers
    plot(x(1),  y(1),  'g^', 'MarkerSize', 11, 'MarkerFaceColor','g', 'DisplayName','Start (t=1)');
    plot(x(60), y(60), 'ks', 'MarkerSize', 11, 'MarkerFaceColor','y', 'DisplayName','Array switch (t=60)');
    plot(x(end),y(end),'rv', 'MarkerSize', 11, 'MarkerFaceColor','r', 'DisplayName','End (t=113)');

    % base station (projected to x-y)
    plot(D.bs(1), D.bs(2), 'kp', 'MarkerSize', 18, 'MarkerFaceColor', [0.5 0.5 0.5], ...
        'DisplayName', sprintf('BS (z=%.1f m)', D.bs(3)));

    axis equal;
    xlim([-halfW-1, halfW+1]); ylim(yl);
    xlabel('x  (lateral, m)'); ylabel('y  (along tunnel, m)');
    title('Tunnel scenario - true vehicle trajectory (top-down)');
    legend('Location','eastoutside');
    set(gca,'FontSize',11);

    if nargin >= 2 && ~isempty(savepath)
        exportgraphics(fig, savepath, 'Resolution', 150);
    end
end
