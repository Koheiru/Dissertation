%% common part for all blocks

% y12 = 1 / 3.5;
% z12 = 0.35;
% k12 = 2.6 * logsig((y12 - 0.1935) * 120.0) - 0.49;
% u01 = 2.6 * logsig(0.1935 * -120);
% u12 = 2.6 * logsig((y12 - 0.1935) * 120.0);
y12 = 0.245;
z12 = 2.6 * 120 * (1 - y12)^2 * logsig((y12 - 0.1935) .* 120.0) * (1.0 - logsig((y12 - 0.1935) .* 120.0));
k12 = 2.6 .* logsig((y12 - 0.1935) .* 120.0) - z12 ./ (1.0 - y12);
u01 = 2.6 * logsig(0.1935 * -120);
u12 = 2.6 .* logsig((y12 - 0.1935) .* 120.0);
f_origin = @(u) (heaviside_restricted(u - u12)) .* (1.0 - z12 ./ (u - k12 + eps)) + ...
                (heaviside_restricted(u - u01) - heaviside_restricted(u - u12)) .* (0.1935 + log(u ./ (2.6 - u) + eps) ./ 120.0);
f_sigm = @(u) logsig(u - 3.0);


%% activation functions
figure();
hold on;
grid on;
xlabel('u');
ylabel('y');

u = -5.0 : 0.01 : 20;
y_sigm = f_sigm(u);
y_origin = f_origin(u);

plot(u, y_sigm,   '--k');
plot(u, y_origin, '-k');
legend('function \sigma', 'function s');
ylim([0 1]);


%% solution curves: activation function variations
figure();
hold on;
grid on;
xlabel('\nu');
ylabel('i');

mu = 0.75;
threshold = 1.0;
theta = 1.0;
alpha = 0.35;%6.5;

f = 0.00 : 0.001 : 0.999;
i = alpha .* f - threshold - mu .* theta .* (log(f) - log(1.0 - f) + 3.0);
plot(f, -i, 'Color', 'black');

f = 0.00 : 0.001 : 0.98;
frequency = f;
middle_frequency = 1.0 / 3.5;
highFilter = heaviside_restricted(frequency - middle_frequency);
lowFilter = heaviside_restricted(frequency) - highFilter;
highActive = 2.46 + 0.35 ./ (1.0 ./ frequency - 1.0);
lowActive  = 2.6 ./ (1.0 + exp((1.0 ./ frequency - 5.2) ./ 0.23));
potential = highFilter .* highActive + lowFilter .* lowActive;
i = alpha .* f - threshold - mu .* theta .* potential;
plot(f, -i, 'Color', 'black');

