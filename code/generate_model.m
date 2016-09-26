function model = generate_model()

  function I = solution_by_input(model)
    I = model.mu .* (model.potential .* model.Q + model.thr - model.frequency .* model.alpha ./ model.mu);
  end

model.mu      = 0.5;
model.mu_min  = 0.00;
model.mu_max  = 1.00;
model.mu_step = 0.01;

model.thr      = 1.00;
model.thr_min  = 0.00;
model.thr_max  = 9.99;
model.thr_step = 0.01;

model.Q      = 1.00;
model.Q_min  = 0.05;
model.Q_max  = 9.00;
model.Q_step = 0.05;

model.I      =  0.00;
model.I_min  = -9.00;
model.I_max  =  9.00;
model.I_step =  0.10;

model.alpha      = 3.0;
model.alpha_min  = 0.00;
model.alpha_max  = 9.99;
model.alpha_step = 0.01;

model.potential = -5.0 : 0.01 : 20.0;
model.frequency = function_f(model.potential);
model.tau = 1.0 ./ model.frequency;

model.activation = @function_f;
model.solution_by_input = @solution_by_input;

model.name = 'Complex';

end

