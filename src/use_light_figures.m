function use_light_figures()
%USE_LIGHT_FIGURES  Force the light MATLAB theme for exported figures.
%
% From R2025a MATLAB figures follow the desktop theme, so a dark-themed
% session (including -batch runs) exports dark-axes PNGs. The report and
% slides expect light figures; call this at the top of every task script.

    try
        s = settings;
        s.matlab.appearance.MATLABTheme.TemporaryValue = 'Light';
    catch
        % pre-R2025a: no theme setting, figures are light by default
    end
end
