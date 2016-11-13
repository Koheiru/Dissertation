%% common part

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
s_origin = @(y) (heaviside_restricted(y - y12)) .* (k12 + z12 ./ (1.0 - y + eps)) + ...
                (heaviside_restricted(y) - heaviside_restricted(y - y12)) .* (2.6 .* logsig(120.0 .* (y - 0.1935)));
dS_origin = @(y) (heaviside_restricted(y - y12)) .* (z12 ./ (1.0 - y + eps) .^ 2) + ...
                 (heaviside_restricted(y) - heaviside_restricted(y - y12)) .* (120.0 * 2.6 .* logsig(120.0 .* (y - 0.1935)) .* (1.0 - logsig(120.0 .* (y - 0.1935))));

f_sigm = @(u) logsig(u - 3.0);
s_sigm = @(y) 3.0 + log(y ./ (1 - y));
dS_sigm = @(y) (1 ./ y + 1.0 ./ (1.0 - y));


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


%%
y12 = 0.245;
z12 = 2.6 * 120 * (1 - y12)^2 * logsig((y12 - 0.1935) .* 120.0) * (1.0 - logsig((y12 - 0.1935) .* 120.0));
k12 = 2.6 .* logsig((y12 - 0.1935) .* 120.0) - z12 ./ (1.0 - y12);
u01 = 2.6 * logsig(0.1935 * -120);
u12 = 2.6 .* logsig((y12 - 0.1935) .* 120.0);
f_origin = @(u) (heaviside_restricted(u - u12)) .* (1.0 - z12 ./ (u - k12 + eps)) + ...
                (heaviside_restricted(u - u01) - heaviside_restricted(u - u12)) .* (0.1935 + log(u ./ (2.6 - u) + eps) ./ 120.0);
s_origin = @(y) (heaviside_restricted(y - y12)) .* (k12 + z12 ./ (1.0 - y + eps)) + ...
                (heaviside_restricted(y) - heaviside_restricted(y - y12)) .* (2.6 .* logsig(120.0 .* (y - 0.1935)));
dS_origin = @(y) (heaviside_restricted(y - y12)) .* (z12 ./ (1.0 - y + eps) .^ 2) + ...
                 (heaviside_restricted(y) - heaviside_restricted(y - y12)) .* (120.0 * 2.6 .* logsig(120.0 .* (y - 0.1935)) .* (1.0 - logsig(120.0 .* (y - 0.1935))));

mu        = 0.75;
threshold = 1.0;
theta     = 1.0;
i         = 2.0;
alpha = 0.01;

y = 0.01 : 0.001 : y12 + 0.1;
dF = alpha  - mu .* theta .* dS_origin(y);

figure();
plot(y, dF); hold on;
grid on;
xlabel('y');
ylabel('dF');


%%
y12 = 0.245;
x_l = exp(-120 * (y12 - 0.1935));
x_r = exp(-120 * (0.0 - 0.1935));
k_l = (x_l + 1)^2 / (156.0 * 2.0 * x_l);
k_r = (x_r + 1)^2 / (156.0 * 2.0 * x_r);

k = 1.0/78.0 : 0.01 : 5.0;
d = (156.0 .* k - 1.0);
x1 = d - sqrt(d.^2 - 1);
x2 = d + sqrt(d.^2 - 1);
y1 = 0.1935 - log(x1) ./ 120.0;
y2 = 0.1935 - log(x2) ./ 120.0;

figure();
plot(k, y1); hold on;
plot(k, y2); hold on;
plot([k_l k_l], [0.0 y12]); hold on;
%plot([k_r k_r], [0.0 y12]); hold on;
plot([min(k) max(k)], [y12 y12]); hold on;
grid on;
xlabel('\mu\theta/\alpha');
ylabel('y');


%%
figure();

y12 = 0.245;
z12 = 2.6 * 120 * (1 - y12)^2 * logsig((y12 - 0.1935) .* 120.0) * (1.0 - logsig((y12 - 0.1935) .* 120.0));
k12 = 2.6 .* logsig((y12 - 0.1935) .* 120.0) - z12 ./ (1.0 - y12);
u01 = 2.6 * logsig(0.1935 * -120);
u12 = 2.6 .* logsig((y12 - 0.1935) .* 120.0);

mu = 0.75;
theta = 1.0;
alpha = 1.0 : 0.01 : 5.0;
y2_high = 1 - sqrt(z12 .* mu .* theta ./ alpha);
i2_high = -alpha .* y2_high + mu .* theta .* (k12 + z12 ./ (1.0 - y2_high));
i1_low = mu .* theta .* u01 .* ones(size(alpha));

subplot(2, 1, 1);
o = mu .* theta ./ alpha;
plot(o, i2_high, '-k'); hold on;
plot(o, i1_low, '-k'); hold on;
grid on;
xlabel('\mu\theta/\alpha');

x = 0.1 : 0.001 : (1 - y12);
o = (x.^2 ./ z12);
y = (k12 - u01) .* x .^2 + 2.0 .* z12 .* x .^ 1 - z12;

k = z12 / (k12 - u01);
o_star = (sqrt(k^2 + k) - k)^2 / z12;

subplot(2, 1, 2);
plot(o, y, '-k'); hold on;
plot([0 0] + o_star, [min(y) max(y)], '--r'); hold on;
grid on;
xlabel('\mu\theta/\alpha');




