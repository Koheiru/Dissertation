function model = generate_model_sigm()

  function frequency = activation(potential)
    frequency = logsig(potential - 3);
  end

  function potential = unactivation(frequency)
    potential = 3 + log(frequency ./ (1.0 - frequency));
  end

  function I = solution_by_input(model)
    I = model.thr - model.alpha .* model.frequency + model.mu .* model.Q .* model.unactivation(model.frequency);
  end

model.mu      = 0.75;
model.mu_min  = 0.00;
model.mu_max  = 1.00;
model.mu_step = 0.01;

model.thr      =  1.00;
model.thr_min  =  0.00;
model.thr_max  = 99.99;
model.thr_step =  0.01;

model.Q      =  1.00;
model.Q_min  =  0.00 + eps;
model.Q_max  = 99.00;
model.Q_step =  0.05;

model.I      =   0.00;
model.I_min  = -99.00;
model.I_max  =  99.00;
model.I_step =   0.10;

model.alpha      =  3.00;
model.alpha_min  =  0.00;
model.alpha_max  = 99.99;
model.alpha_step =  0.01;

model.potential = -5.0 : 0.01 : 20.0;
model.frequency = activation(model.potential);
model.tau = 1.0 ./ model.frequency;

model.activation = @activation;
model.unactivation = @unactivation;
model.solution_by_input = @solution_by_input;

model.name = 'Sigmoid';

end

