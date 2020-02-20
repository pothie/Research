function h = paper_figure(varargin) % (1)
% similar to the MATLAB function figure , it can use the same input arguments
hsize = 8.5; vsize = 11; border = .4; % (2)
paper_position = [0 0 hsize vsize] + border*[1 1 -2 -2]; % (3)
ratio = hsize/vsize; % (4)
screen_size = get(0, 'ScreenSize'); % (5)
window_size_p(4) = (2/3)*screen_size(4); % (6)
window_size_p(3) = ratio*window_size_p(4); % (7)
window_size_p(2) = 1; % (8)
window_size_p(1) = screen_size(3) - window_size_p(3); % (9)
all_h = get(groot, 'Children'); % (10)
h = figure(varargin{1:end}); % (11)
if ~ismember(h, all_h) % (12)
    set(h, 'Position', window_size_p, ...
    'PaperPositionMode', 'manual', ...
    'PaperPosition', paper_position); % (13)
end % (14)
if nargout == 0 % (15)
clear h % (16)
end % (17)
end