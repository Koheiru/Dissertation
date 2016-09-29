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

f_sigm = @(u) logsig(u - 3.0);
s_sigm = @(y) 3.0 + log(y ./ (1 - y));


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

% y_star = zeros(1, length(pre_alpha));
% for n = 1 : length(pre_alpha)
%   y_star(n) = fsolve(@(y) pre_alpha(n) * y - pre_i(n) - threshold - mu * theta * s_sigm(y), 0.5);
% end
% 
% y1_star = zeros(1, length(post_alpha));
% y2_star = zeros(1, length(post_alpha));
% for n = 1 : length(post_alpha)
%   y1_star(n) = fsolve(@(y) post_alpha(n) * y - post_i(n) - threshold - mu * theta * s_sigm(y), 0.1);
%   y2_star(n) = fsolve(@(y) post_alpha(n) * y - post_i(n) - threshold - mu * theta * s_sigm(y), 0.9);
% end

figure_subplot(1, 2, 1);
hold on; grid off; box on;
xlabel('\alpha');
ylabel('y^{*}');
plot(pre_alpha, y_star, '.k');
plot(post_alpha, y1_star, '.k');
plot(post_alpha, y2_star, '.k');
ylim([-0.05 1.05]);

figure_subplot(1, 2, 2);
hold on; grid off; box on;
xlabel('i');
ylabel('\alpha');
plot(-i1, alpha, '-k');
plot(-i2, alpha, '-k');
plot([-pre_i, -post_i], [pre_alpha post_alpha], ':k');

%% [sigmoidal model] bifurcation: F surface
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


%% [sigmoidal model] solution curves: alpha variations
figure();
hold on;
grid on;
xlabel('\nu');
ylabel('i');


alpha = 2.5;
threshold = 1.0;
theta = 1.0;
for mu = 0.1 : 0.1 : 0.9
  f = 0.01 : 0.01 : 0.99;
  i = alpha .* f - threshold - mu .* theta .* (log(f) - log(1.0 - f) + 3.0);
  plot(f, -i, 'Color', 'black');
end

% mu = 0.75;
% threshold = 1.0;
% theta = 1.0;
% for alpha = [0 1 2.5 5 10 15]
%   f = 0.01 : 0.01 : 0.99;
%   i = alpha .* f - threshold - mu .* theta .* (log(f) - log(1.0 - f) + 3.0);
%   plot(f, -i, 'Color', 'black');
% end


%% [original model]: replace s- and f-function derivative!
% %% [original model] bifurcation: equlibrium intervals
% figure();
% hold on;
% grid on;
% 
% x_bound = 1.0 / 3.5;
% x = 0.01 : 0.01 : x_bound;
% k = 1.0 ./ (1.0 + exp((1.0 ./ x - 5.2) ./ 0.23));
% y = 2.6 .* k .* (1.0 - k) ./ (0.23 .* x.^2);
% 
% plot(x, y);
% 
% x = x_bound : 0.01 : 0.95;
% y = 0.35 ./ (x - 1.0).^2;
% 
% plot(x, y);
% 
% mu = 0.75;
% threshold = 1.0;
% theta = 1.0;
% alpha = 6.5;
% 
% plot([0.0 1.0], [alpha./(mu.*theta) alpha./(mu.*theta)]);
% 
% 
% %% [original model] bifurcation: i vs alpha diagram
% figure();
% hold on;
% grid on;
% xlabel('i');
% ylabel('\alpha');
% 
% mu = 0.75;
% threshold = 1.0;
% theta = 1.0;
% 
% %alpha = 1.4 * mu * theta / 3.0 : 0.1 : 25;
% alpha = 1.4 * mu * theta / 3.0  .* [1 2 3];
% f1 = (1.0 + sqrt(-1.4 .* mu .* theta ./ alpha + 3.0)) ./ 2.0;
% i1 = alpha .* f1 - threshold - mu .* theta .* (2.11 + 0.35 ./ (1.0 - f1));
% f2 = (1.0 - sqrt(-1.4 .* mu .* theta ./ alpha + 3.0)) ./ 2.0;
% i2 = alpha .* f2 - threshold - mu .* theta .* (2.11 + 0.35 ./ (1.0 - f2));
% 
% plot(i1, alpha, 'Color', 'green');
% plot(i2, alpha, 'Color', 'blue');










