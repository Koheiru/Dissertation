%
% Activation function with non-linear saturation.
%
function frequency = function_f(potential)

y12 = 0.245;
z12 = 2.6 * 120 * (1 - y12)^2 * logsig((y12 - 0.1935) .* 120.0) * (1.0 - logsig((y12 - 0.1935) .* 120.0));
k12 = 2.6 .* logsig((y12 - 0.1935) .* 120.0) - z12 ./ (1.0 - y12);
u01 = 2.6 * logsig(0.1935 * -120);
u12 = 2.6 .* logsig((y12 - 0.1935) .* 120.0);

high_filter = heaviside_restricted(potential - u12);
low_filter  = heaviside_restricted(potential - u01) - high_filter;

high_active = 1.0 - z12 ./ (potential - k12 + eps);
low_active  = 0.1935 + log(potential ./ (2.6 - potential) + eps) ./ 120.0;

frequency  = high_filter .* high_active + low_filter .* low_active;
