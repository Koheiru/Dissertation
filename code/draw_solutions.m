function draw_solutions()

  function drawSolution(model)
    cla;
    hold on;

    title(['Ensemble equilibrium (', model.name, ')']);
    xlabel('\nu');
    ylabel('I');
    grid on;

    L = model.solution_by_input(model);
    plot(model.frequency, L, 'Color', 'black', 'LineWidth', 2);

    G_values = zeros(size(L)) + model.I;
    plot(model.frequency, G_values, 'Color', 'red', 'LineWidth', 2);

    ylim([-3.0 4.0]);
    xlim([-0.01 0.99]);
  end
  
  function drawInputDependency(model)
    cla;
    hold on;

    title(['Input dependency (', model.name, ')']);
    xlabel('I');
    ylabel('\nu');
    grid on;
    
    I = -0.1 : 0.01 : 6.0;
    L = model.solution_by_input(model);

    L1 = zeros(length(I), 1);
    L2 = zeros(length(I), 1);
    for i = 1 : length(I)
      indexes = find(abs(L - I(i)) < 0.01);
      if (~isempty(indexes))
        L1(i) = model.frequency(indexes(1));
        L2(i) = model.frequency(indexes(end));
      end
    end

    plot(I, L1, 'Color', 'black', 'LineWidth', 2);
    plot(I, L2, '--', 'Color', 'black', 'LineWidth', 2);

    plot([model.I model.I], [0 1], 'Color', 'red', 'LineWidth', 2);
  end
      
  function drawProperties(handle)
    figure(handle);
    models = get(handle, 'UserData');
    
    width = 5;
    height = length(models) + 1;
    
    for index = 1 : length(models)
      model = models{index};
      offset = index * width;
      subplot(height, width, offset + 2 : offset + 3); drawSolution(model);
      subplot(height, width, offset + 4 : offset + 5); drawInputDependency(model);
    end
  end

  function onISpinChanged(spinner, event)
    models = get(gcf, 'UserData');
    
    for index = 1 : length(models)
      models{index}.I = spinner.getValue();
    end
    
    set(gcf, 'UserData', models);
    drawProperties(gcf);
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

      drawProperties(gcf);
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
    
    width = 5;
    height = length(models) + 1;
    
    subplot(height, width, 1);
    axis off;
    
    panel = uipanel('Position', get(gca, 'Position'));
    container = uigridcontainer('v0', 'Parent', panel, 'Units', 'norm', 'Position', [.1,.1,.9,.9]);
    set(container, 'GridSize', [1,2]);
    
    label_obj = javax.swing.JLabel('I =');
    javacomponent(label_obj, [0 0 1 1], container);
    param_model = javax.swing.SpinnerNumberModel(0.0, -5, 25, 0.1);

    spin_obj = javax.swing.JSpinner(param_model);
    spin_obj.setName('I');
    spin_handle = javacomponent(spin_obj, [0 0 1 1], container);
    set(spin_handle, 'MouseWheelMovedCallback', @onSpinScrolled);
    set(spin_handle, 'StateChangedCallback', @onISpinChanged);
        
    for index = 1 : length(models)
      offset = index * width;
      
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
%models_{3} = generate_model_lin();

figure_handle_ = 500 + randi(5000);
prepareFigure(figure_handle_, models_);
drawProperties(figure_handle_);

end

