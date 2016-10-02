function handle = figure_subplot(h, w, p)
  if (length(p) > 1)
    error('Fix me, friend!');
  end
  
  ncols = w;
  nrows = h;
  row = floor((p - 1) / ncols) + 1;
  col = mod(p - 1, ncols) + 1;
  
  if (ncols == 3)
    axisw = (1 / ncols) * 0.74;
    axisl = 0.075 + (axisw+0.075) * (col - 1);
  elseif (ncols == 2)
    axisw = 0.35;
    axisl = 0.1 + 0.5 * (col - 1);
	elseif (ncols == 1)
    axisw = 0.85;
    axisl = 0.1;
  end
  
  if (nrows == 4)
    axish = (1 / nrows) * 0.70;
    axisb = 0.025 + (axish+0.065) * (row - 1);
  elseif (nrows == 2)
    axish = (1 / nrows) * 0.70;
    axisb = 0.05 + (axish+0.125) * (row - 1);
  elseif (nrows == 1)
    axish = 0.75;
    axisb = 0.075;
  end

  handle = subplot(h, w, p, 'position', [axisl, 1.0 - axisb - axish, axisw, axish]);
  
  set(handle, 'FontName', 'Times New Roman');
  set(handle, 'FontSize', 10);

  xhandle = get(handle,'XLabel');
  set(xhandle, 'FontName', 'Times New Roman');
  set(xhandle, 'FontSize', 12);

  yhandle = get(handle,'YLabel');
  set(yhandle, 'FontName', 'Times New Roman');
  set(yhandle, 'FontAngle', 'italic');
  set(yhandle, 'FontSize', 12);
  set(yhandle, 'Rotation', 0.0);
end

