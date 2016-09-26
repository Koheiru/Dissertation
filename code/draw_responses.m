function draw_responses()

  function f = activationByDistance(model, R)
    I = 8.0 .* exp(-R ./ 10.0);
    f = model.activation((I - model.thr) / model.Q);
  end

  function drawResponseSurface(model, w, x, y, z)
    cla;
    hold on;

    title(['Response surface (', model.name, ')']);
    xlabel('x_{1}');
    ylabel('x_{2}');
    zlabel('\nu');
    grid on;
    
    surf(x, y, z);
    
    xlim([min(x) max(x)]);
    ylim([min(y) max(y)]);
    zlim([0.0 1.0]);
    view(-40, 20);
  end

  function drawResponseCountour(model, w, x, y, z)
    cla;
    hold on;

    title(['Response surface (', model.name, ')']);
    xlabel('x_{1}');
    ylabel('x_{2}');
    grid on;
    
    contour(x, y, z, 0.1 : 0.1 : 0.9);
    plot(w(1), w(2), 'r*');
  end
      
  function drawResponses(handle)
    figure(handle);
    models = get(handle, 'UserData');
    
    width = 7;
    height = length(models);
    
    target = [0.0 0.0];
    x1_values = -30 : 0.5 : 30;
    x2_values = -30 : 0.5 : 30;
    
    points = zeros(length(x1_values) * length(x2_values), 2);
    count = length(x1_values);
    for i = 1 : count
      tmp = [ones(length(x1_values), 1) * x1_values(i), x2_values'];
      points((i - 1) * count + 1 : i * count, :) = tmp(:, :);
    end
    
    distance = pdist2(target, points, 'euclidean');
    
    for index = 1 : length(models)
      model = models{index};
      offset = (index - 1) * width;
      
      activation = activationByDistance(model, distance);
      activation = reshape(activation, length(x1_values), length(x2_values));
      
      subplot(height, width, offset + 2 : offset + 4); drawResponseSurface(model, target, x1_values, x2_values, activation);
      subplot(height, width, offset + 5 : offset + 7); drawResponseCountour(model, target, x1_values, x2_values, activation);
    end
  end

  function onSpinChanged(index, spinner, event)
    models = get(gcf, 'UserData');
    operator = java.lang.StringBuilder().append('models{index}.').append(spinner.getName());
    value_getter = operator.toString();

    value = spinner.getValue();
    if value ~= eval(value_getter)
      operator.append(' = value;');
      value_setter = operator.toString();

      eval(value_setter);
      set(gcf, 'UserData', models);

      drawResponses(gcf);
    end
  end

  function onSpinScrolled(spinner, event)
    if event.getWheelRotation() > 0
      value = spinner.getModel().getPreviousValue();
    else
      value = spinner.getModel().getNextValue();
    end
    spinner.setValue(value);
  end

  function appendParamSpinner(container, param_name, index)
    models = get(gcf, 'UserData');
    
    label_obj = javax.swing.JLabel([param_name, ' =']);
    javacomponent(label_obj, [0 0 1 1], container);

    param_value = eval(['models{index}.', param_name]);
    param_min   = eval(['models{index}.', param_name, '_min']);
    param_max   = eval(['models{index}.', param_name, '_max']);
    param_step  = eval(['models{index}.', param_name, '_step']);
    param_model = javax.swing.SpinnerNumberModel(param_value, param_min, param_max, param_step);

    spin_obj = javax.swing.JSpinner(param_model);
    spin_obj.setName(param_name);
    spin_handle = javacomponent(spin_obj, [0 0 1 1], container);
    set(spin_handle, 'MouseWheelMovedCallback', @onSpinScrolled);
    set(spin_handle, 'StateChangedCallback', @(h,e)onSpinChanged(index,h,e));
  end

  function prepareFigure(figure_handle, models)
    figure(figure_handle);
    set(figure_handle, 'UserData', models);
    
    width = 7;
    height = length(models);
    
    for index = 1 : length(models)
      offset = (index - 1) * width;
      
      subplot(height, width, offset + 1);
      axis off;
      
      panel = uipanel('Position', get(gca, 'Position'));
      container = uigridcontainer('v0', 'Parent', panel, 'Units', 'norm', 'Position', [.1,.1,.9,.9]);

      set(container, 'GridSize', [4,2]);
      appendParamSpinner(container, 'Q', index);
      appendParamSpinner(container, 'mu', index);
      appendParamSpinner(container, 'alpha', index);
      appendParamSpinner(container, 'thr', index);
    end
  end


models_ = cell(1, 2);
models_{1} = generate_model();
models_{2} = generate_model_sigm();

figure_handle_ = 654;
prepareFigure(figure_handle_, models_);
drawResponses(figure_handle_);

end

