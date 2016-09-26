%
% Invertion of activation function (see f_function).
%
function potential = function_s(frequency)

y12 = 0.245;
z12 = 2.6 * 120 * (1 - y12)^2 * logsig((y12 - 0.1935) .* 120.0) * (1.0 - logsig((y12 - 0.1935) .* 120.0));
k12 = 2.6 .* logsig((y12 - 0.1935) .* 120.0) - z12 ./ (1.0 - y12);
u01 = 2.6 * logsig(0.1935 * -120);
u12 = 2.6 .* logsig((y12 - 0.1935) .* 120.0);

high_filter = heaviside_restricted(frequency - y12);
low_filter = heaviside_restricted(frequency) - high_filter;

high_active = k12 + z12 ./ (1.0 - frequency + eps);
low_active  = 2.6 .* logsig(120.0 .* (frequency - 0.1935));

potential = high_filter .* high_active + low_filter .* low_active;
