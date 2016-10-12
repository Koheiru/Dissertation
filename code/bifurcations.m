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
s_origin = @(y) (heaviside_restricted(y - y12)) .* (k12 + z12 ./ (1.0 - y + eps)) + ...
                (heaviside_restricted(y) - heaviside_restricted(y - y12)) .* (2.6 .* logsig(120.0 .* (y - 0.1935)));
dS_origin = @(y) (heaviside_restricted(y - y12)) .* (z12 ./ (1.0 - y + eps) .^ 2) + ...
                 (heaviside_restricted(y) - heaviside_restricted(y - y12)) .* (120.0 * 2.6 .* logsig(120.0 .* (y - 0.1935)) .* (1.0 - logsig(120.0 .* (y - 0.1935))));

f_sigm = @(u) logsig(u - 3.0);
s_sigm = @(y) 3.0 + log(y ./ (1 - y));
dS_sigm = @(y) (1 ./ y + 1.0 ./ (1.0 - y));


%% [sigmoidal model] bifurcation: i vs alpha diagram
fh = figure();
figure_adjust(fh, [17.5 6.5]);

mu = 0.75;
threshold = 1.0;
theta = 1.0;
alpha = 4.0 * mu * theta : 0.1 : 20;

f1 = (1.0 + sqrt(1.0 - 4.0 .* mu .* theta ./ alpha)) ./ 2.0;
f2 = (1.0 - sqrt(1.0 - 4.0 .* mu .* theta ./ alpha)) ./ 2.0;
i1 = alpha .* f1 - threshold - mu .* theta .* (log(f1) - log(1.0 - f1) + 3.0);
i2 = alpha .* f2 - threshold - mu .* theta .* (log(f2) - log(1.0 - f2) + 3.0);

i_mid = (i1 + i2) ./ 2.0;
k = (i_mid(2) - i_mid(1)) / (alpha(2) - alpha(1));
pre_alpha = 0 : 0.1 : 4.0 * mu * theta;
pre_i = pre_alpha .* k + (i_mid(1) - pre_alpha(end) .* k);
pre_alpha = pre_alpha(1 : end - 1);
pre_i = pre_i(1 : end - 1);
post_alpha = alpha;
post_i = i_mid;

y_star = zeros(1, length(pre_alpha));
for n = 1 : length(pre_alpha)
  y0 = 0.1; y_star(n) = find_solution_dyn(pre_alpha(n), -pre_i(n), threshold, mu, theta, f_sigm, y0, theta * s_sigm(y0));
end

y1_star = zeros(1, length(post_alpha));
y2_star = zeros(1, length(post_alpha));
for n = 1 : length(post_alpha)
  y0 = 0.1; y1_star(n) = find_solution_dyn(post_alpha(n), -post_i(n), threshold, mu, theta, f_sigm, y0, theta * s_sigm(y0));
  y0 = 0.9; y2_star(n) = find_solution_dyn(post_alpha(n), -post_i(n), threshold, mu, theta, f_sigm, y0, theta * s_sigm(y0));
end

figure_subplot(1, 2, 1);
hold on; grid off; box on;
xlabel('\alpha');
ylabel('y^{*}');
plot(pre_alpha, y_star, '.k');
plot(post_alpha, y1_star, '.k');
plot(post_alpha, y2_star, '.k');
ylim([-0.05 1.05]);

figure_subplot(1, 2, 2);
plot(alpha, -i1, '-k'); hold on;
plot(alpha, -i2, '-k'); hold on;
plot([pre_alpha post_alpha], [-pre_i, -post_i], ':k'); hold on;
grid off; box on;
xlabel('\alpha');
ylabel('i');

%% [sigmoidal model] bifurcation: i vs theta diagram
fh = figure();
figure_adjust(fh, [17.5 6.5]);

mu = 0.75;
threshold = 1.0;
alpha = 4.5;
theta = [0.01 : 0.05 : alpha / (4.0 * mu) - eps, alpha / (4.0 * mu)];

f1 = (1.0 + sqrt(1.0 - 4.0 .* mu .* theta ./ alpha)) ./ 2.0;
f2 = (1.0 - sqrt(1.0 - 4.0 .* mu .* theta ./ alpha)) ./ 2.0;
i1 = alpha .* f1 - threshold - mu .* theta .* (log(f1) - log(1.0 - f1) + 3.0);
i2 = alpha .* f2 - threshold - mu .* theta .* (log(f2) - log(1.0 - f2) + 3.0);

i_mid = (i1 + i2) ./ 2.0;
k = (i_mid(end) - i_mid(end-1)) / (theta(end) - theta(end-1));
post_theta = alpha / (4.0 * mu) : 0.1 : 3.0;
post_i = post_theta .* k + (i_mid(end) - post_theta(1) .* k);
pre_theta = theta(1 : end - 1);
pre_i = i_mid(1 : end - 1);

y_star = zeros(1, length(post_theta));
for n = 1 : length(post_theta)
  y0 = 0.5; y_star(n) = find_solution_dyn(alpha, -post_i(n), threshold, mu, post_theta(n), f_sigm, y0, post_theta(n) * s_sigm(y0));
end

y1_star = zeros(1, length(pre_theta));
y2_star = zeros(1, length(pre_theta));
for n = 1 : length(pre_theta)
  y0 = 0.1; y1_star(n) = find_solution_dyn(alpha, -pre_i(n), threshold, mu, pre_theta(n), f_sigm, y0, pre_theta(n) * s_sigm(y0));
  y0 = 0.9; y2_star(n) = find_solution_dyn(alpha, -pre_i(n), threshold, mu, pre_theta(n), f_sigm, y0, pre_theta(n) * s_sigm(y0));
end

figure_subplot(1, 2, 1);
plot(post_theta, y_star, '.k'); hold on;
plot(pre_theta, y1_star, '.k'); hold on;
plot(pre_theta, y2_star, '.k'); hold on;
grid off; box on;
xlabel('\theta');
ylabel('y^{*}');
ylim([-0.05 1.05]);

figure_subplot(1, 2, 2);
plot(theta, -i1, '-k'); hold on;
plot(theta, -i2, '-k'); hold on;
plot([pre_theta, post_theta], [-pre_i, -post_i], ':k'); hold on;
grid off; box on;
xlabel('\theta');
ylabel('i');

%% [sigmoidal model] bifurcation: F-surface defined by alpha
fh = figure();
figure_adjust(fh, [17.5 10.5]);

mu = 0.75;
threshold = 1.0;
theta = 1.0;
alpha = 2.0 : 0.5 : 15;
y = [0.00001, 0.0001, 0.001, 0.01:0.01:0.1, 0.1 : 0.1 : 0.9, 0.9:0.01:0.99, 0.999, 0.9999, 0.99999];
i = -(alpha' * y - threshold - repmat(mu .* theta .* s_sigm(y), length(alpha), 1));

alpha_s = 4.0 * mu * theta : 0.5 : 15;
f1 = (1.0 + sqrt(1.0 - 4.0 .* mu .* theta ./ alpha_s)) ./ 2.0;
f2 = (1.0 - sqrt(1.0 - 4.0 .* mu .* theta ./ alpha_s)) ./ 2.0;
i1 = alpha_s .* f1 - threshold - mu .* theta .* (log(f1) - log(1.0 - f1) + 3.0);
i2 = alpha_s .* f2 - threshold - mu .* theta .* (log(f2) - log(1.0 - f2) + 3.0);

surf(y, alpha, i, 'FaceColor', 'white');
grid on; hold on;
plot3(f1, alpha_s, -i1, '-k', 'LineWidth', 2);
plot3(f2, alpha_s, -i2, '-k', 'LineWidth', 2);
xlabel('y');
ylabel('\alpha');
zlabel('i');
view([-141, 36]);

%% [sigmoidal model] bifurcation: F-surface defined by theta
fh = figure();
figure_adjust(fh, [17.5 10.5]);

mu = 0.75;
threshold = 1.0;
alpha = 4.5;
theta = 0.5 : 0.1 : 2.0;
y = [0.00001, 0.0001, 0.001, 0.01:0.01:0.1, 0.1 : 0.1 : 0.9, 0.9:0.01:0.99, 0.999, 0.9999, 0.99999];
i = -(repmat(alpha .* y, length(theta), 1) - threshold - mu .* theta' * s_sigm(y));

theta_s = [theta(1) : (theta(end) - theta(1)) / length(theta) : alpha / (4.0 * mu), alpha / (4.0 * mu)];
f1 = (1.0 + sqrt(1.0 - 4.0 .* mu .* theta_s ./ alpha)) ./ 2.0;
f2 = (1.0 - sqrt(1.0 - 4.0 .* mu .* theta_s ./ alpha)) ./ 2.0;
i1 = alpha .* f1 - threshold - mu .* theta_s .* (log(f1) - log(1.0 - f1) + 3.0);
i2 = alpha .* f2 - threshold - mu .* theta_s .* (log(f2) - log(1.0 - f2) + 3.0);

surf(y, theta, i, 'FaceColor', 'white');
grid on; hold on;
plot3(f1, theta_s, -i1, '-k', 'LineWidth', 2);
plot3(f2, theta_s, -i2, '-k', 'LineWidth', 2);
xlabel('y');
ylabel('\theta');
zlabel('i');
view([-154, 28]);

%% [sigmoidal model] bifurcation: i-surfaces defined by alpha and theta
fh = figure();
figure_adjust(fh, [17.5 8.5]);

mu = 0.75;
threshold = 1.0;
theta = 0.1 : 0.1 : 3.5;
alpha = 0.1 : 0.5 : 10;

f1 = (1.0 + sqrt(1.0 - 4.0 .* mu .* theta' * (1.0 ./ alpha))) ./ 2.0;
f2 = (1.0 - sqrt(1.0 - 4.0 .* mu .* theta' * (1.0 ./ alpha))) ./ 2.0;
f1((1.0 ./ theta') * alpha < 4.0 .* mu) = NaN; f1 = real(f1);
f2((1.0 ./ theta') * alpha < 4.0 .* mu) = NaN; f2 = real(f2);
i1 = repmat(alpha, length(theta), 1) .* f1 - threshold - mu .* repmat(theta', 1, length(alpha)) .* (log(f1) - log(1.0 - f1) + 3.0);
i2 = repmat(alpha, length(theta), 1) .* f2 - threshold - mu .* repmat(theta', 1, length(alpha)) .* (log(f2) - log(1.0 - f2) + 3.0);

surf(alpha, theta, -i1, 'FaceColor', 'w'); hold on;
surf(alpha, theta, -i2, 'FaceColor', 'w'); hold on;
plot3(alpha, alpha ./ (4.0 * mu), -(alpha .* 0.5 - threshold - mu .* (alpha ./ (4.0 * mu)) .* 3.0), '-k', 'LineWidth', 2); hold on;
grid on; box off;
xlabel('\alpha');
ylabel('\theta');
zlabel('i');
view([117 33]);



%% [original model] bifurcation: i vs alpha diagram
fh = figure();
figure_adjust(fh, [17.5 6.5]);

mu = 0.75;
threshold = 1.0;
theta = 1.0;
alpha = [0.0 : 0.5 : 65];

x_l = exp(-120 * (y12 - 0.1935));
x_r = exp(-120 * (0.0 - 0.1935));
k_l = (x_l + 1)^2 / (156.0 * 2.0 * x_l);
k_r = (x_r + 1)^2 / (156.0 * 2.0 * x_r);
d = (156.0 .* mu .* theta ./ alpha - 1);

f2_low = 0.0 + + zeros(size(alpha));
i2_low = mu * theta * u01 - threshold + zeros(size(alpha));

x1_low = d + sqrt(d.^2 - 1);
x1_low(78.0 .* mu .* theta ./ alpha < 1.0) = NaN;
x1_low(mu .* theta ./ alpha > (x_r + 1.0)^2/(312.0*x_r)) = NaN; 
x1_low = real(x1_low);
f1_low = 0.1935 - log(x1_low) / 120.0;
i1_low = alpha .* f1_low - threshold - mu .* theta .* (2.6 .* logsig(120.0 .* (f1_low - 0.1935)));

f1_high = 1.0 - sqrt(z12 .* mu .* theta ./ alpha);
f1_high(z12 .* mu .* theta ./ alpha > (1.0 - y12)^2) = NaN; 
f1_high = real(f1_high);
i1_high = alpha .* f1_high - threshold - mu .* theta .* (k12 + z12 ./ (1.0 - f1_high));

x2_high = d - sqrt(d.^2 - 1);
x2_high(78.0 .* mu .* theta ./ alpha < 1.0) = NaN; 
x2_high(mu .* theta ./ alpha > (x_l + 1.0)^2/(312.0*x_l)) = NaN; 
x2_high = real(x2_high);
f2_high = 0.1935 - log(x2_high) / 120.0;
i2_high = alpha .* f2_high - threshold - mu .* theta .* (2.6 .* logsig(120.0 .* (f2_high - 0.1935)));

plot(alpha, -i2_low, '-g'); hold on;
plot(alpha, -i1_low, '-g'); hold on;
plot(alpha, -i2_high, '-r'); hold on;
plot(alpha, -i1_high, '-r'); hold on;
grid off; box on;
xlabel('\alpha');
ylabel('i');
set(fh, 'Name', ['theta = ', num2str(theta)]);

%% [original model] bifurcation: i vs theta diagram
fh = figure();
figure_adjust(fh, [17.5 6.5]);

mu = 0.75;
threshold = 1.0;
theta = [0.0 : 0.001 : 5.0];
alpha = 5.0;

x_l = exp(-120 * (y12 - 0.1935));
x_r = exp(-120 * (0.0 - 0.1935));
k_l = (x_l + 1)^2 / (156.0 * 2.0 * x_l);
k_r = (x_r + 1)^2 / (156.0 * 2.0 * x_r);
d = (156.0 .* mu .* theta ./ alpha - 1);

f2_low = 0.0 + zeros(size(theta));
i2_low = mu * theta * u01 - threshold + zeros(size(theta));

x1_low = d + sqrt(d.^2 - 1);
x1_low(78.0 .* mu .* theta ./ alpha < 1.0) = NaN;
x1_low(mu .* theta ./ alpha > (x_r + 1.0)^2/(312.0*x_r)) = NaN; 
x1_low = real(x1_low);
f1_low = 0.1935 - log(x1_low) / 120.0;
i1_low = alpha .* f1_low - threshold - mu .* theta .* (2.6 .* logsig(120.0 .* (f1_low - 0.1935)));

f1_high = 1.0 - sqrt(z12 .* mu .* theta ./ alpha);
f1_high(z12 .* mu .* theta ./ alpha > (1.0 - y12)^2) = NaN; 
f1_high = real(f1_high);
i1_high = alpha .* f1_high - threshold - mu .* theta .* (k12 + z12 ./ (1.0 - f1_high));

x2_high = d - sqrt(d.^2 - 1);
x2_high(78.0 .* mu .* theta ./ alpha < 1.0) = NaN; 
x2_high(mu .* theta ./ alpha > (x_l + 1.0)^2/(312.0*x_l)) = NaN; 
x2_high = real(x2_high);
f2_high = 0.1935 - log(x2_high) / 120.0;
i2_high = alpha .* f2_high - threshold - mu .* theta .* (2.6 .* logsig(120.0 .* (f2_high - 0.1935)));

plot(theta, -i2_low, '-g'); hold on;
plot(theta, -i1_low, '-g'); hold on;
plot(theta, -i2_high, '-r'); hold on;
plot(theta, -i1_high, '-r'); hold on;
grid off; box on;
xlabel('\theta');
ylabel('i');
set(fh, 'Name', ['alpha = ', num2str(alpha)]);

%% [original model] bifurcation: F surface defined by alpha
fh = figure();
figure_adjust(fh, [17.5 10.5]);

mu = 0.75;
threshold = 1.0;
theta = 1.0;
alpha = [0.0 : 0.1 : 1.0, 1.0 : 0.5 : 30];
y = 0.01 : 0.01 : 0.98;
i = -(alpha' * y - threshold - repmat(mu .* theta .* s_origin(y), length(alpha), 1));

% copy-paste from above: i vs alpha diagram
x_l = exp(-120 * (y12 - 0.1935));
x_r = exp(-120 * (0.0 - 0.1935));
k_l = (x_l + 1)^2 / (156.0 * 2.0 * x_l);
k_r = (x_r + 1)^2 / (156.0 * 2.0 * x_r);
d = (156.0 .* mu .* theta ./ alpha - 1);

f2_low = 0.0 + + zeros(size(alpha));
i2_low = mu * theta * u01 - threshold + zeros(size(alpha));

x1_low = d + sqrt(d.^2 - 1);
x1_low(78.0 .* mu .* theta ./ alpha < 1.0) = NaN;
x1_low(mu .* theta ./ alpha > (x_r + 1.0)^2/(312.0*x_r)) = NaN; 
x1_low = real(x1_low);
f1_low = 0.1935 - log(x1_low) / 120.0;
i1_low = alpha .* f1_low - threshold - mu .* theta .* (2.6 .* logsig(120.0 .* (f1_low - 0.1935)));

f1_high = 1.0 - sqrt(z12 .* mu .* theta ./ alpha);
f1_high(z12 .* mu .* theta ./ alpha > (1.0 - y12)^2) = NaN; 
f1_high = real(f1_high);
i1_high = alpha .* f1_high - threshold - mu .* theta .* (k12 + z12 ./ (1.0 - f1_high));

x2_high = d - sqrt(d.^2 - 1);
x2_high(78.0 .* mu .* theta ./ alpha < 1.0) = NaN; 
x2_high(mu .* theta ./ alpha > (x_l + 1.0)^2/(312.0*x_l)) = NaN; 
x2_high = real(x2_high);
f2_high = 0.1935 - log(x2_high) / 120.0;
i2_high = alpha .* f2_high - threshold - mu .* theta .* (2.6 .* logsig(120.0 .* (f2_high - 0.1935)));

surf(y, alpha, i, 'FaceColor', 'white'); hold on; grid on;
plot3(f1_low, alpha, -i1_low, '-g', 'LineWidth', 2);
plot3(f2_low, alpha, -i2_low, '-g', 'LineWidth', 2);
plot3(f1_high, alpha, -i1_high, '-r', 'LineWidth', 2);
plot3(f2_high, alpha, -i2_high, '-r', 'LineWidth', 2);
xlabel('y');
ylabel('\alpha');
zlabel('i');
view([-141, 36]);

%% [original model] bifurcation: F surface defined by theta
fh = figure();
figure_adjust(fh, [17.5 10.5]);

mu = 0.75;
threshold = 1.0;
alpha = 2.5;
theta = 0.0 : 0.1 : 6.0;
y = 0.01 : 0.01 : 0.98;
i = -(repmat(alpha .* y, length(theta), 1) - threshold - mu .* theta' * s_origin(y));

% copy-paste from above: i vs theta diagram
x_l = exp(-120 * (y12 - 0.1935));
x_r = exp(-120 * (0.0 - 0.1935));
k_l = (x_l + 1)^2 / (156.0 * 2.0 * x_l);
k_r = (x_r + 1)^2 / (156.0 * 2.0 * x_r);
d = (156.0 .* mu .* theta ./ alpha - 1);

f2_low = 0.0 + zeros(size(theta));
i2_low = mu * theta * u01 - threshold + zeros(size(theta));

x1_low = d + sqrt(d.^2 - 1);
x1_low(78.0 .* mu .* theta ./ alpha < 1.0) = NaN;
x1_low(mu .* theta ./ alpha > (x_r + 1.0)^2/(312.0*x_r)) = NaN; 
x1_low = real(x1_low);
f1_low = 0.1935 - log(x1_low) / 120.0;
i1_low = alpha .* f1_low - threshold - mu .* theta .* (2.6 .* logsig(120.0 .* (f1_low - 0.1935)));

f1_high = 1.0 - sqrt(z12 .* mu .* theta ./ alpha);
f1_high(z12 .* mu .* theta ./ alpha > (1.0 - y12)^2) = NaN; 
f1_high = real(f1_high);
i1_high = alpha .* f1_high - threshold - mu .* theta .* (k12 + z12 ./ (1.0 - f1_high));

x2_high = d - sqrt(d.^2 - 1);
x2_high(78.0 .* mu .* theta ./ alpha < 1.0) = NaN; 
x2_high(mu .* theta ./ alpha > (x_l + 1.0)^2/(312.0*x_l)) = NaN; 
x2_high = real(x2_high);
f2_high = 0.1935 - log(x2_high) / 120.0;
i2_high = alpha .* f2_high - threshold - mu .* theta .* (2.6 .* logsig(120.0 .* (f2_high - 0.1935)));

surf(y, theta, i, 'FaceColor', 'white');
grid on; hold on;
plot3(f1_low, theta, -i1_low, '-g', 'LineWidth', 2);
plot3(f2_low, theta, -i2_low, '-g', 'LineWidth', 2);
plot3(f1_high, theta, -i1_high, '-r', 'LineWidth', 2);
plot3(f2_high, theta, -i2_high, '-r', 'LineWidth', 2);
xlabel('y');
ylabel('\theta');
zlabel('i');
view([-154, 28]);

