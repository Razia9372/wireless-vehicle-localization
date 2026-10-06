function outfile = make_animation(D, est, label, savepath, opts)
%MAKE_ANIMATION  Top-down video of true vs estimated position with rays.
%
%   outfile = make_animation(D, est, label, savepath, opts)
%
% Draws, per time step: the tunnel, BS, the true position, the estimated
% position with a fading trail, and the active multipath rays (measured
% range + measured angles of the active array). Writes an MP4 (falls back
% to Motion-JPEG AVI if the MPEG-4 profile is unavailable).
%
% opts: .fps (default 15), .trail (default 20 steps)
% Returns the actual output filename written.

    if nargin < 5, opts = struct; end
    fps   = getf(opts,'fps',15);
    trail = getf(opts,'trail',20);

    [vw, outfile] = open_writer(savepath, fps);
    open(vw); cleaner = onCleanup(@() close(vw));

    fig = figure('Color','w','Position',[100 100 1100 420],'Visible','off');
    ax = axes(fig); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
    halfW = max(abs(D.ue_pos(1,:)))+0.8;
    xlim(ax,[min(D.ue_pos(2,:))-4, max(D.ue_pos(2,:))+4]); ylim(ax,[-halfW halfW]);
    xlabel(ax,'y  (along tunnel, m)'); ylabel(ax,'x  (lateral, m)');
    set(ax,'FontSize',11);

    for t = 1:D.N
        cla(ax);
        % tunnel walls + BS
        plot(ax, xlim(ax), [ halfW-0.3  halfW-0.3],'k-','LineWidth',1.2);
        plot(ax, xlim(ax), [-halfW+0.3 -halfW+0.3],'k-','LineWidth',1.2);
        plot(ax, D.bs(2), 0, 'kp','MarkerSize',16,'MarkerFaceColor',[.5 .5 .5]);

        % active rays (measured), array-local -> global
        k = D.ant_idx(t); az0 = D.ant_or(k,1); el0 = D.ant_or(k,2);
        hRay = gobjects(1,0);
        for i = 1:numel(D.toa_hat{t})
            r = D.toa_hat{t}(i);
            gAz = D.aoa_hat{t}(1,i)+az0; gEl = D.aoa_hat{t}(2,i)+el0;
            ep = D.bs + r*[cosd(gEl)*cosd(gAz); cosd(gEl)*sind(gAz); sind(gEl)];
            h = plot(ax, [D.bs(2) ep(2)], [0 ep(1)], '-', 'Color',[1 .7 .2 .6], 'LineWidth',1);
            if i == 1, hRay = h; end
        end

        % trail + current estimate + truth
        i0 = max(1,t-trail);
        hTrail = plot(ax, est(2,i0:t), est(1,i0:t), '-', 'Color',[0 .45 .74], 'LineWidth',1.6);
        hTrue  = plot(ax, D.ue_pos(2,t), D.ue_pos(1,t), 'ko','MarkerSize',9,'MarkerFaceColor','k');
        hEst   = plot(ax, est(2,t), est(1,t), 'o','MarkerSize',9,'MarkerFaceColor',[0 .45 .74],'MarkerEdgeColor','w');

        title(ax, sprintf('%s   t=%d/%d   array %d   err=%.2f m', label, t, D.N, k, ...
            hypot(est(1,t)-D.ue_pos(1,t), est(2,t)-D.ue_pos(2,t))));
        legend(ax, [hRay hTrail hTrue hEst], ...
            [repmat({'measured rays'},1,numel(hRay)), {'est trail','true','estimate'}], ...
            'AutoUpdate','off','Location','eastoutside');
        drawnow;
        writeVideo(vw, getframe(fig));
    end
    close(fig);
end

function [vw, outfile] = open_writer(savepath, fps)
    [p,n,~] = fileparts(savepath);
    try
        outfile = fullfile(p,[n '.mp4']);
        vw = VideoWriter(outfile,'MPEG-4');
    catch
        outfile = fullfile(p,[n '.avi']);
        vw = VideoWriter(outfile,'Motion JPEG AVI');
    end
    vw.FrameRate = fps;
end

function v = getf(s,f,d)
    if isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
