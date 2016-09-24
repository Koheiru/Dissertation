function figure_export(figure_handle, target_size)

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

set(0,'defaultAxesFontName', 'Times New Roman');
set(0,'defaultTextFontName', 'Times New Roman');

set(0,'units','centimeters');
screen_scale = get(0,'screensize') ./ [1 1 target_width target_height];
set(0,'units','pixels');
screen_size = get(0,'screensize');
frame_position = zeros(1, 4);
frame_position(3) = screen_size(3) / screen_scale(3);
frame_position(4) = screen_size(4) / screen_scale(4);
frame_position(1) = screen_size(3) - frame_position(3);
frame_position(2) = 0.0;

fh = figure();
set(fh, 'Units', 'pixels');
set(fh, 'Position', frame_position);
%showplottool('propertyeditor');

objects = allchild(target);
copyobj(objects, fh);
%close(figure_handle);

set(fh, 'PaperUnits', 'centimeters');
set(fh, 'PaperSize', [target_width target_height]);
set(fh, 'PaperPosition', [0.0 0.0 target_width target_height]);

ah_list = findobj(fh, 'type', 'axes');
ah_count = length(ah_list);

% horizontal subplots
%{
if (ah_count > 2)
  error('Fix me, friend!');
elseif (ah_count > 1)
  set(ah_list(1), 'Units', 'normalized');
  set(ah_list(1), 'Position', [0.1 0.175 0.35 0.8]);
  set(ah_list(2), 'Units', 'normalized');
  set(ah_list(2), 'Position', [0.6 0.175 0.35 0.8]);
elseif (ah_count > 0)
  set(ah_list(1), 'Units', 'normalized');
  set(ah_list(1), 'Position', [0.1 0.175 0.85 0.8]);
end
%}

% vertical subplots
%{
if (ah_count > 3)
  error('Fix me, friend!');
elseif (ah_count > 2)
  set(ah_list(1), 'Units', 'normalized');
  set(ah_list(2), 'Units', 'normalized');
  set(ah_list(3), 'Units', 'normalized');
  set(ah_list(1), 'Position', [0.1 (0.075 + 0.1)         0.85 0.175]);
  set(ah_list(2), 'Position', [0.1 (0.35 + 0.025 + 0.1)  0.85 0.175]);
  set(ah_list(3), 'Position', [0.1 (0.65  + 0.025 + 0.1)  0.85 0.175]);
elseif (ah_count > 1)
  set(ah_list(1), 'Units', 'normalized');
  set(ah_list(2), 'Units', 'normalized');
  set(ah_list(1), 'Position', [0.1 (0.075 + 0.1)        0.85 0.3]);
  set(ah_list(2), 'Position', [0.1 (0.5 + 0.025 + 0.1)  0.85 0.3]);
elseif (ah_count > 0)
  set(ah_list(1), 'Units', 'normalized');
  set(ah_list(1), 'Position', [0.1 0.175 0.85 0.8]);
end
%}

% squared subplots


for ah = ah_list(:)'
  set(ah, 'FontName', 'Times New Roman');
  set(ah, 'FontSize', 10);

  xh = get(ah,'XLabel');
  set(xh, 'FontName', 'Times New Roman');
  set(xh, 'FontSize', 12);

  yh = get(ah,'YLabel');
  set(yh, 'FontName', 'Times New Roman');
  set(yh, 'FontAngle', 'italic');
  set(yh, 'FontSize', 12);
  set(yh, 'Rotation', 0.0);
end

filename = ['figures/figure_', num2str(randi(1000))];
print('-deps', '-r300', '-tiff', '-painters', filename);

close(fh);
