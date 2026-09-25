function save_png(f, file, resolution)
%SAVE_PNG  Export a figure to PNG robustly.
%
%   save_png(f, file, resolution)
%
%   Writes to a temporary file first and then moves it into place, with a
%   few retries. Overwriting a PNG that another program (image viewer,
%   Explorer preview, indexer) has open otherwise fails intermittently with
%   "PNG library failed". The axes hover toolbar is hidden in the image.

if nargin < 3, resolution = 200; end
axs = findall(f, 'Type', 'axes');
for a = axs.', a.Toolbar.Visible = 'off'; end
cleanup = onCleanup(@() restore(axs));

tmp = [tempname '.png'];
lastErr = [];
for attempt = 1:4
    try
        exportgraphics(f, tmp, 'Resolution', resolution);
        [ok, msg] = movefile(tmp, file, 'f');
        if ok, return; end
        lastErr = MException('save_png:Move', '%s', msg);
    catch err
        lastErr = err;
    end
    pause(0.5*attempt);
end
if isfile(tmp), delete(tmp); end
warning('save_png:Failed', 'Could not save %s (%s). Close any program that has it open and run again.', ...
    file, lastErr.message);
end

function restore(axs)
for a = axs.'
    if isvalid(a), a.Toolbar.Visible = 'on'; end
end
end
