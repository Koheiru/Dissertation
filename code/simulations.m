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


%% [sigmoidal model] simulations in various regimes
fh = figure();
figure_adjust(fh, [26.2 14.0]);

duration  = 1000;

u         = zeros(duration, 2);
y         = zeros(duration, 2);
theta     = zeros(duration, 1) + 1.0;
i         = zeros(duration, 1);
mu        = [0.75   0.75];
threshold = [1.0    1.0];
alpha     = [1.5    3.5];

% 1st demo
i( 26 : 100) = 2.0 .* [1 : 75] ./ 75;
i(100 : 126) = 2.0;
i(126 : 200) = 2.0 .* (1.0 - [1 : 75] ./ 75);

% 2nd demo
    i(226 : 400) = 1.0;
theta(226 : 300) = 1.0 - 0.25 .* (0.0 + [1 : 75] ./ 75);
theta(300 : 326) = 0.75;
theta(326 : 400) = 1.0 - 0.25 .* (1.0 - [1 : 75] ./ 75);

decision = 100;
dt = 1.0 / decision;
for t = 1 : duration - 1
  u(t + 1, :) = u(t, :) + dt .* (alpha .* y(t, :) + i(t) - threshold - mu .* u(t, :));
  for k = 2 : decision
    u(t + 1, :) = u(t + 1, :) + dt .* (alpha .* y(t, :) + i(t) - threshold - mu .* u(t + 1, :));
  end
  y(t + 1, :) = [f_sigm(u(t + 1, 1) / theta(t)), f_sigm(u(t + 1, 2) / theta(t))];
end

f_plus  = (1.0 + sqrt(1.0 - 4.0 .* mu(2) .* theta ./ alpha(2))) ./ 2.0;
f_minus = (1.0 - sqrt(1.0 - 4.0 .* mu(2) .* theta ./ alpha(2))) ./ 2.0;
i_plus  = -alpha(2) .* f_plus  + threshold(2) + mu(2) .* theta .* s_sigm(f_plus);
i_minus = -alpha(2) .* f_minus + threshold(2) + mu(2) .* theta .* s_sigm(f_minus);

subplot(3, 1, 1);
plot([1 : duration], y(:, 1), 'r-'); hold on;
plot([1 : duration], y(:, 2), 'b--'); hold on;
grid off;
xlabel('t');
ylabel('y');
ylim([-0.05 1.05]);

subplot(3, 1, 2);
plot([1 : duration], i, 'k-'); hold on;
plot([1 : duration], i_plus, 'k--'); hold on;
plot([1 : duration], i_minus, 'k--'); hold on;
grid off;
xlabel('t');
ylabel('i');

subplot(3, 1, 3);
plot([1 : duration], theta, 'k-'); hold on;
grid off;
xlabel('t');
ylabel('\theta');
















