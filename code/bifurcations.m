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
figure();
hold on;
grid on;
xlabel('i');
ylabel('\alpha');

mu = 0.75;
threshold = 1.0;
theta = 1.0;
alpha = 4.0 * mu * theta : 0.1 : 25;

f1 = (1.0 + sqrt(1.0 - 4.0 .* mu .* theta ./ alpha)) ./ 2.0;
i1 = alpha .* f1 - threshold - mu .* theta .* (log(f1) - log(1.0 - f1) + 3.0);

f2 = (1.0 - sqrt(1.0 - 4.0 .* mu .* theta ./ alpha)) ./ 2.0;
i2 = alpha .* f2 - threshold - mu .* theta .* (log(f2) - log(1.0 - f2) + 3.0);

plot(-i1, alpha, 'Color', 'green');
plot(-i2, alpha, 'Color', 'blue');


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










