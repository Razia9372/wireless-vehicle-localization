function fig = plot_traj_overlay(D, ests, names, ttl, savepath)
%PLOT_TRAJ_OVERLAY  Overlay estimated trajectories on the true one (top-down).
%
%   fig = plot_traj_overlay(D, ests, names, ttl, savepath)
%
% Inputs:
%   ests  - cell array of 2xN estimates [x;y]
%   names - cell array of legend names (same length as ests)
%   ttl   - title string
%   savepath - optional PNG path
%
% The tunnel is long (y~100 m) and narrow (x~+/-3 m); axes are NOT equal so
% the lateral behaviour is visible. True path is thick black; estimates are
% coloured. BS and array-switch are marked.

    x = D.ue_pos(1,:); y = D.ue_pos(2,:);
    co = lines(max(numel(ests),3));

    fig = figure('Color','w','Position',[100 100 1000 520]);
    hold on; box on; grid on;

    % tunnel walls (style of the brief's example result figure)
    halfW = max(abs(x)) + 0.5;
    plot([min(y)-5, max(y)+5], [ halfW  halfW], 'k-', 'LineWidth', 1.5, 'HandleVisibility','off');
    plot([min(y)-5, max(y)+5], [-halfW -halfW], 'k-', 'LineWidth', 1.5, 'HandleVisibility','off');

    plot(y, x, 'k-', 'LineWidth', 2.4, 'DisplayName', 'true');
    for i = 1:numel(ests)
        plot(ests{i}(2,:), ests{i}(1,:), '-', 'Color', co(i,:), ...
            'LineWidth', 1.4, 'DisplayName', names{i});
    end
    plot(D.bs(2), 0, 'kp', 'MarkerSize', 16, 'MarkerFaceColor',[.5 .5 .5], ...
        'DisplayName','BS (x=0)');
    xline(D.ue_pos(2,60), 'k--', 'array switch', 'HandleVisibility','off', ...
        'LabelVerticalAlignment','bottom');

    xlabel('y  (along tunnel, m)'); ylabel('x  (lateral, m)');
    ylim([-6 6]);
    title(ttl); legend('Location','northoutside','Orientation','horizontal');
    set(gca,'FontSize',11);

    if nargin >= 5 && ~isempty(savepath)
        exportgraphics(fig, savepath, 'Resolution', 150);
    end
end
