function figure_adjust(figure_handle, target_size)

target = gcf;
if (nargin > 0)
  target = figure_handle;
end

% A4 full is 21.0x29.7, bordered is 17.5x25.7
target_width  = 17.5;
target_height = 24.0;
if (nargin > 1)
  target_width  = target_size(1);
  target_height = target_size(2);
end

set(0,'units','centimeters');
screen_scale = get(0,'screensize') ./ [1 1 target_width target_height];
set(0,'units','pixels');
screen_size = get(0,'screensize');
frame_position = get(target, 'Position');
frame_position(2) = frame_position(2) - (screen_size(4) / screen_scale(4) - frame_position(4));
frame_position(3) = screen_size(3) / screen_scale(3);
frame_position(4) = screen_size(4) / screen_scale(4);

set(target, 'Units', 'pixels');
set(target, 'Position', frame_position);
