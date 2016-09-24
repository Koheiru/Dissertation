%% common part for all blocks

% k12 = 2.6 * logsig(38.73 / 3.5) - 0.49;
% u01 = 2.6 * logsig(-23.22);
% u12 = 2.6 * logsig(38.73 / 3.5);
% f_origin = @(u) heaviside_restricted(u - u12) .* (1.0 - 0.35 ./ (u - k12)) + ...
%                 (heaviside_restricted(u - u01) - heaviside_restricted(u - u12)) .* (0.1935 + log(u ./ (2.6 - u)) ./ 120.0);

f_sigm = @(u) logsig(u - 3.0);
s_sigm = @(y) 3.0 + log(y ./ (1 - y));

%% [sigmoidal model] figures of equilibrium: soft and hard states
figure();

mu = 0.75;
threshold = 1.25;
i_soft = 3.0;
i_hard = 1.0;
alpha_soft = 1.0;
alpha_hard = 5.0;
theta = 1.0;

y = 0.001 : 0.001 : 0.999;
F_soft = i_soft + alpha_soft .* y - threshold - mu .* theta .* s_sigm(y);
I_soft = -(alpha_soft .* y - threshold - mu .* theta .* s_sigm(y));
F_hard = i_hard + alpha_hard .* y - threshold - mu .* theta .* s_sigm(y);
I_hard = -(alpha_hard .* y -threshold - mu .* theta .* s_sigm(y));

subplot(2, 2, 1);
hold on; grid off; box on;
plot(I_soft, y, '-k');
plot([i_soft i_soft], [0 1], '--k');
xlabel('i');
ylabel('y');
xlim([-1 5]);
ylim([-0.02 1.02]);

subplot(2, 2, 3);
hold on; grid off; box on;
plot(y, F_soft, '-k');
plot([0 1], [0 0], '--k');
xlabel('y');
ylabel('F(y)');
ylim([-3.0 3.0])

subplot(2, 2, 2);
hold on; grid off; box on;
plot(I_hard, y, '-k');
plot([i_hard i_hard], [0 1], '--k');
xlabel('i');
ylabel('y');
xlim([-1 4]);
ylim([-0.02 1.02]);

subplot(2, 2, 4);
hold on; grid off; box on;
plot(y, F_hard, '-k');
plot([0 1], [0 0], '--k');
xlabel('y');
ylabel('F(y)');
ylim([-3.0 3.0])


%% [sigmoidal model] figures of equilibrium conditions
figure();

mu = 0.75;
threshold = 1.0;
theta     = 1.0;
i         = 1.0;
alpha_zero = 4.0 * mu * theta;
alpha_soft = alpha_zero / 1.5;
alpha_hard = alpha_zero * 1.5;

y = 0.001 : 0.001 : 0.999;
F_zero = i + alpha_zero .* y - threshold - mu .* theta .* s_sigm(y);
F_soft = i + alpha_soft .* y - threshold - mu .* theta .* s_sigm(y);
F_hard = i + alpha_hard .* y - threshold - mu .* theta .* s_sigm(y);

subplot(2, 3, 1);
hold on; grid off; box on;
plot(y, F_soft, '-k');
plot([0 1], [0 0], '--k');
xlabel('y');
ylabel('F(y)');

subplot(2, 3, 2);
hold on; grid off; box on;
plot(y, F_zero, '-k');
plot([0 1], [0 0], '--k');
xlabel('y');
ylabel('F(y)');

subplot(2, 3, 3);
hold on; grid off; box on;
plot(y, F_hard, '-k');
plot([0 1], [0 0], '--k');
xlabel('y');
ylabel('F(y)');

y = 0.1 : 0.01 : 0.9;
dF_zero = alpha_zero - mu .* theta .* (1 ./ y + 1.0 ./ (1.0 - y));
dF_soft = alpha_soft - mu .* theta .* (1 ./ y + 1.0 ./ (1.0 - y));
dF_hard = alpha_hard - mu .* theta .* (1 ./ y + 1.0 ./ (1.0 - y));

subplot(2, 3, 4);
hold on; grid off; box on;
plot(y, dF_soft, '-k');
plot([0 1], [0 0], '--k');
xlabel('y');
ylabel('F''(y)');
ylim([-3 3]);

subplot(2, 3, 5);
hold on; grid off; box on;
plot(y, dF_zero, '-k');
plot([0 1], [0 0], '--k');
xlabel('y');
ylabel('F''(y)');
ylim([-3 3]);

subplot(2, 3, 6);
hold on; grid off; box on;
plot(y, dF_hard, '-k');
plot([0 1], [0 0], '--k');
xlabel('y');
ylabel('F''(y)');
ylim([-3 3]);


%%
figure();

mu = 0.75;
threshold = 1.0;
alpha = 2.5;
theta = 1.0;

subplot(2, 2, 1);
hold on; grid on;
xlabel('y');
ylabel('\alpha');
zlabel('i');

y = 0.001 : 0.01 : 0.999;
alpha_values = 0.0 : 0.1 : 10.0;
i = -alpha_values' * y + repmat(mu .* theta .* s_sigm(y), length(alpha_values), 1);
surf(y, alpha_values, i, 'FaceColor', 'none');
view(-132, 24);

%%
figure();

mu = 0.75;
threshold = 1.0;
alpha = 2.5;
theta = 1.0;

hold on; grid on;
xlabel('y');
ylabel('\theta');
zlabel('i');

mu = 0.75;
threshold = 1.0;
alpha = 4.5;
theta = 1.0;

y = 0.001 : 0.01 : 0.999;
theta_values = 0.01 : 0.1 : 2.0;
i = repmat(-alpha .* y, length(theta_values), 1) + mu .* theta_values' * s_sigm(y);
surf(y, theta_values, i, 'FaceColor', 'none');
view(-132, 24);










