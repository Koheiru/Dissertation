function [y, u] = find_solution_dyn(alpha, i, threshold, mu, theta, f, y0, u0)

  params.dt             = 0.01;
  params.max_iterations = 200 / params.dt;
  params.max_delta      = 0.000000001;

dt              = params.dt;
max_iterations  = params.max_iterations;
max_delta       = params.max_delta;

u_cur = u0;
y_cur = y0;
for k = 2 : max_iterations
  u_prev = u_cur;
  y_prev = y_cur;
  
  du = alpha * y_prev + i - threshold - mu * u_prev;
  u_cur = u_prev + dt * du;
  y_cur = f(u_cur / theta);
  
  delta  = abs(y_cur - y_prev);
  if (delta < max_delta)
    break;
  end
end

y = y_cur;
u = u_cur;

end

