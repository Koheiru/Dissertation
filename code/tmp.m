%%
figure();
hold on; grid on;

mu = 0.75;
threshold = 1.0;
theta = 1.0;
alpha = 10.0; % 58.5
input = 0.0;
y = 0.000 : 0.001 : 0.99;

frequency = y;
middle_frequency = 1.0 / 3.5;
high_filter = heaviside_restricted(frequency - middle_frequency);
low_filter = heaviside_restricted(frequency) - high_filter;
high_active = 2.11 + 0.35 ./ (1.0 - frequency);
low_active  = 2.6 .* logsig((frequency - 0.1935) .* 120.0);
potential = high_filter .* high_active + low_filter .* low_active;
g = potential;

F = alpha .* y + input - threshold - mu .* theta .* g;
plot(y, -F, 'Color', 'black', 'LineWidth', 1);

%%
%y1 = 1.0 - sqrt(0.35 .* mu .* theta ./ alpha)
d = 1.0 - 1.3 .* 120.0 .* mu  .* theta ./ alpha;
y2 = 0.1935 - log(-d - sqrt(d .^ 2 - 1.0)) ./ 120.0
y3 = 0.1935 - log(-d + sqrt(d .^ 2 - 1.0)) ./ 120.0


%%
figure();
hold on;
grid on;
xlabel('\nu');
ylabel('u');

x = 0.01 : 0.001 : 1.0 - 0.01;
frequency = x;
tau = 1.0 ./ (frequency);
middle_frequency = 1.0 / 3.5;
high_filter = heaviside_restricted(frequency - middle_frequency);
low_filter = heaviside_restricted(frequency) - high_filter;
high_active = 2.46 + 0.35 ./ (tau - 1.0);
low_active  = 2.6 ./ (1.0 + exp((tau - 5.2) ./ 0.23));
potential = high_filter .* high_active + low_filter .* low_active;
y = potential;
plot(x, y, 'Color', 'blue', 'LineWidth', 1);

x = 0.01 : 0.001 : 1.0 - 0.01;
frequency = x;
middle_frequency = 1.0 / 3.5;
high_filter = heaviside_restricted(frequency - middle_frequency);
low_filter = heaviside_restricted(frequency) - high_filter;
high_active = 2.11 + 0.35 ./ (1.0 - frequency);
low_active  = 2.6 .* logsig((frequency - 0.1935) .* 120.0);
potential = high_filter .* high_active + low_filter .* low_active;
y = potential;
plot(x, y, 'Color', 'blue', 'LineWidth', 1);



%%
f = 0;
y1 = 2.6 .* logsig((f - 0.1935) .* 120.0);
y2 = 2.11 + 0.35 ./ (1.0 - f)






%%
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
k12 = 2.6 * logsig(38.73/3.5) - 0.49;
u01 = 2.6 * logsig(-23.22);
u12 = 2.6 * logsig(38.73/3.5);
highFilter = heaviside_restricted(potential - 2.6);
lowFilter  = heaviside_restricted(potential - 0.0) - highFilter;
highActive = 1.0 - 0.35 ./ (potential - k12);
lowActive  = 0.1935 + log(potential ./ (2.6 - potential)) ./ 120.0;
frequency  = highFilter .* highActive + lowFilter .* lowActive;
y2 = frequency;

subplot(3, 2, 1);
hold on; grid on;
xlabel('u');
ylabel('y');
plot(x, y1, 'Color', 'black', 'LineWidth', 1);
ylim([0 1]);

subplot(3, 2, 3);
hold on; grid on;
xlabel('u');
ylabel('y');
plot(x, y2, 'Color', 'blue', 'LineWidth', 1);
ylim([0 1]);

subplot(3, 2, 5);
hold on; grid on;
xlabel('u');
ylabel('\Delta');
plot(x, abs(y1 - y2), 'Color', 'black', 'LineWidth', 1);



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
middle_frequency = 1.0 / 3.5;
high_filter = heaviside_restricted(frequency - middle_frequency);
low_filter = heaviside_restricted(frequency) - high_filter;
high_active = 2.6 .* logsig(38.73/3.5) - 0.49 + 0.35 ./ (1.0 - frequency);
low_active  = 2.6 .* logsig((frequency - 0.1935) .* 120.0);
potential = high_filter .* high_active + low_filter .* low_active;
y2 = potential;

subplot(3, 2, 2);
hold on; grid on;
xlabel('y');
ylabel('u');
plot(x, y1, 'Color', 'black', 'LineWidth', 1);
ylim([0 10]);

subplot(3, 2, 4);
hold on; grid on;
xlabel('y');
ylabel('u');
plot(x, y2, 'Color', 'blue', 'LineWidth', 1);
ylim([0 10]);

subplot(3, 2, 6);
hold on; grid on;
xlabel('y');
ylabel('\Delta u');
plot(x, abs(y1 - y2), 'Color', 'black', 'LineWidth', 1);


%%
figure();

mu = 0.75;
alpha = 60.0;
theta = 0.0 : 0.1 : alpha/(4*mu);

f1 = 0.5 - sqrt(0.25 - mu .* theta ./ alpha);
f2 = 0.5 + sqrt(0.25 - mu .* theta ./ alpha);
i1 = -alpha .* f1 + mu .* theta .* (3.0 + log(f1 ./ (1.0 - f1)));
i2 = -alpha .* f2 + mu .* theta .* (3.0 + log(f2 ./ (1.0 - f2)));

z1 = sqrt(0.25 - mu .* theta ./ alpha);
K = @(n,x) (2 ^ (2 * n + 3)) .* (x .^ (2 * n + 3)) ./ (4 * n^2 + 8 * n + 3);
ii_tail = K(0,z1) + K(1,z1) + K(2,z1) + K(3,z1) + K(4,z1) + K(5,z1);
ii1 = alpha .* (0.25 - ii_tail);
ii2 = alpha .* (0.25 + ii_tail);

z2 = 2.0 .* sqrt(0.25 - mu .* theta ./ alpha);
K = @(n,x) (x .^ (2 * n + 3)) ./ (4 * n^2 + 8 * n + 3);
iii_tail = K(0,z2) + K(1,z2) + K(2,z2) + K(3,z2) + K(4,z2) + K(5,z2);
iii1 = 0.25 .* alpha .* (1.0 - 4.* iii_tail);
iii2 = 0.25 .* alpha .* (1.0 + 4.* iii_tail);

subplot(2, 1, 1);
plot(z1, ii1, '-b'); hold on;
plot(z1, ii2, '-b'); hold on;
grid on;

subplot(2, 1, 2);
plot(theta, i1, '-r', 'LineWidth', 4); hold on;
plot(theta, i2, '-r', 'LineWidth', 4); hold on;
plot(theta, ii1, '-b', 'LineWidth', 4); hold on;
plot(theta, ii2, '-b', 'LineWidth', 4); hold on;
plot(theta, iii1, '-g'); hold on;
plot(theta, iii2, '-g'); hold on;
grid on;


%%
figure();

A = 0.1;
K = @(n,x) (log(A / (1.0 - A)) / (A * (1.0 - A)^n)) .* (x - A).^n;
x = 0.0 : 0.01 : 1.0;

y0 = log(x ./ (1.0 - x));
plot(x, y0, '-r'); hold on;

for n = 2000
  y = K(0,x);
  for k = 1 : n
    y = y + K(k,x);
  end
  %plot(x, abs(y - y0), '-b'); hold on;
  plot(x, -y, '-b'); hold on;
end






