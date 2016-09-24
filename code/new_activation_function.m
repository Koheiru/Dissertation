%% smooth s-function: parameters lookup
figure();

skew = 120.0;
x = 0.01 : 0.001 : 0.3;
y1 = 2.6 ./ (1.0 + exp((1.0 ./ x - 5.2) ./ 0.23));
y2 = 2.6 .* logsig((x - 0.1935) .* skew);

subplot(2, 1, 1);
hold on; grid on;
plot(x, y1, 'Color', 'black', 'LineWidth', 1);
plot(x, y2, 'Color', 'blue', 'LineWidth', 1);

subplot(2, 1, 2);
hold on; grid on;
title(['skew is ', num2str(skew)]);
ylabel('\Delta y');
plot(x, abs(y1 - y2), 'Color', 'black', 'LineWidth', 1);

%% smooth s-function: sub-functions plot
figure();

x = 0.01 : 0.001 : 1.0 - 0.01;
frequency = x;
% y12 = 1.0 / 3.5;
% z12 = 0.35;
% k12 = 2.6 .* logsig(38.73 / 3.5) - 0.49;
y12 = 0.245;
z12 = 2.6 * 120 * (1 - y12)^2 * logsig((y12 - 0.1935) .* 120.0) * (1.0 - logsig((y12 - 0.1935) .* 120.0));
k12 = 2.6 .* logsig((y12 - 0.1935) .* 120.0) - z12 ./ (1.0 - y12);
high_filter = heaviside_restricted(frequency - y12);
low_filter = heaviside_restricted(frequency) - high_filter;
high_active = k12 + z12 ./ (1.0 - frequency + eps);
low_active  = 2.6 .* logsig((frequency - 0.1935) .* 120.0);
potential = high_filter .* high_active + low_filter .* low_active;
y = potential;

subplot(1, 2, 1);
hold on; grid on;
xlabel('y');
ylabel('u');
plot(x, low_active, 'Color', 'black', 'LineWidth', 1);
plot(x, high_active, 'Color', 'black', 'LineWidth', 1);
ylim([0 10]);

subplot(1, 2, 2);
hold on; grid on;
xlabel('y');
ylabel('u');
plot(x, y, 'Color', 'black', 'LineWidth', 1);
ylim([0 10]);


%% f-function vs smooth f-function
figure();

x = -2 : 0.001 : 10;

potential = x;
highFilter = heaviside_restricted(potential - 2.6 + eps);
lowFilter  = heaviside_restricted(potential) - highFilter;
highActive = 1.0 + 0.35 ./ (potential - 2.46) + eps;
lowActive  = 5.2 + 0.23 .* log(2.6 ./ potential - 1.0) + eps;
frequency  = highFilter ./ highActive + lowFilter ./ lowActive;
y1 = frequency;

potential = x;
y12 = 0.245;
z12 = 2.6 * 120 * (1 - y12)^2 * logsig((y12 - 0.1935) .* 120.0) * (1.0 - logsig((y12 - 0.1935) .* 120.0));
k12 = 2.6 .* logsig((y12 - 0.1935) .* 120.0) - z12 ./ (1.0 - y12);
u01 = 2.6 * logsig(0.1935 * -120);
u12 = 2.6 .* logsig((y12 - 0.1935) .* 120.0);
% y12 = 1 / 3.5;
% z12 = 0.35;
% k12 = 2.6 * logsig((y12 - 0.1935) * 120.0) - 0.49;
% u01 = 2.6 * logsig(0.1935 * -120);
% u12 = 2.6 * logsig((y12 - 0.1935) * 120.0);
highFilter = heaviside_restricted(potential - u12);
lowFilter  = heaviside_restricted(potential - u01) - highFilter;
highActive = 1.0 - z12 ./ (potential - k12 + eps);
lowActive  = 0.1935 + log(potential ./ (2.6 - potential) + eps) ./ 120.0;
frequency  = highFilter .* highActive + lowFilter .* lowActive;
y2 = frequency;

subplot(3, 1, 1);
hold on; grid on;
xlabel('u');
ylabel('old');
plot(x, y1, 'Color', 'black', 'LineWidth', 1);
ylim([0 1]);

subplot(3, 1, 2);
hold on; grid on;
xlabel('u');
ylabel('new');
plot(x, y2, 'Color', 'black', 'LineWidth', 1);
ylim([0 1]);

subplot(3, 1, 3);
hold on; grid on;
xlabel('u');
ylabel('\Delta');
plot(x, abs(y1 - y2), 'Color', 'black', 'LineWidth', 1);


%% s-function vs smooth s-function
figure();

x = 0.01 : 0.001 : 1.0 - 0.01;

frequency = x;
tau = 1.0 ./ (frequency);
middle_frequency = 1.0 / 3.5;
high_filter = heaviside_restricted(frequency - middle_frequency);
low_filter = heaviside_restricted(frequency) - high_filter;
high_active = 2.46 + 0.35 ./ (tau - 1.0);
low_active  = 2.6 ./ (1.0 + exp((tau - 5.2) ./ 0.23));
potential = high_filter .* high_active + low_filter .* low_active;
y1 = potential;

frequency = x;
% y12 = 1.0 / 3.5;
% z12 = 0.35;
% k12 = 2.6 .* logsig(38.73 / 3.5) - 0.49;
y12 = 0.245;
z12 = 2.6 * 120 * (1 - y12)^2 * logsig((y12 - 0.1935) .* 120.0) * (1.0 - logsig((y12 - 0.1935) .* 120.0));
k12 = 2.6 .* logsig((y12 - 0.1935) .* 120.0) - z12 ./ (1.0 - y12);
high_filter = heaviside_restricted(frequency - y12);
low_filter = heaviside_restricted(frequency) - high_filter;
high_active = k12 + z12 ./ (1.0 - frequency + eps);
low_active  = 2.6 .* logsig((frequency - 0.1935) .* 120.0);
potential = high_filter .* high_active + low_filter .* low_active;
y2 = potential;

subplot(3, 1, 1);
hold on; grid on;
xlabel('y');
ylabel('u');
plot(x, y1, 'Color', 'black', 'LineWidth', 1);
ylim([0 10]);

subplot(3, 1, 2);
hold on; grid on;
xlabel('y');
ylabel('u');
plot(x, y2, 'Color', 'blue', 'LineWidth', 1);
ylim([0 10]);

subplot(3, 1, 3);
hold on; grid on;
xlabel('y');
ylabel('\Delta u');
plot(x, abs(y1 - y2), 'Color', 'black', 'LineWidth', 1);
ylim([0 0.2]);





