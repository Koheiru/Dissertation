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

%% [sigmoidal model] differential models: model simulation
mu = 0.75;
threshold = 1.0;
alpha = 3.5;
theta = 1.0;

f = @(u) f_sigm(u);
du = @(u, y, i) alpha .* y + i - threshold - mu .* u;

duration = 175;
i = zeros(1, duration);
i( 25 :  50) = 0.5;
i( 50 :  75) = 1.2;
i( 75 : 100) = 3.5;
i(100 : 125) = 0.5;
i(125 : 150) = 0.0;

dt_e10 = 10;
u_e10 = zeros(1, duration * dt_e10);
y_e10 = zeros(1, duration * dt_e10);

dt_e100 = 100;
u_e100 = zeros(1, duration * dt_e100);
y_e100 = zeros(1, duration * dt_e100);

dt_rk10 = 10;
u_rk10 = zeros(1, duration * dt_rk10);
y_rk10 = zeros(1, duration * dt_rk10);

dt_pk10 = 10;
u_pk10 = zeros(1, duration * dt_pk10);
y_pk10 = zeros(1, duration * dt_pk10);

for t = 2 : duration
  for k = 1 : dt_e10
    T = (t - 1) * dt_e10;
    u_e10(T + k) = u_e10(T + k - 1) + du(u_e10(T + k - 1), y_e10(T + k - 1), i(t)) / dt_e10;
    y_e10(T + k) = f(u_e10(T + k) / theta);
  end
  
  for k = 1 : dt_e100
    T = (t - 1) * dt_e100;
    u_e100(T + k) = u_e100(T + k - 1) + du(u_e100(T + k - 1), y_e100(T + k - 1), i(t)) / dt_e100;
    y_e100(T + k) = f(u_e100(T + k) / theta);
  end
  
  for k = 1 : dt_rk10
    T = (t - 1) * dt_rk10;
    h = 1.0 / dt_rk10;
    k1 = du(u_rk10(T + k - 1), y_rk10(T + k - 1), i(t));
    k2 = du(u_rk10(T + k - 1) + k1 * h / 2.0, f((u_rk10(T + k - 1) + k1 * h / 2.0) / theta), i(t));
    k3 = du(u_rk10(T + k - 1) + k2 * h / 2.0, f((u_rk10(T + k - 1) + k2 * h / 2.0) / theta), i(t));
    k4 = du(u_rk10(T + k - 1) + k3 * h, f((u_rk10(T + k - 1) + k3 * h) / theta), i(t));
    u_rk10(T + k) = u_rk10(T + k - 1) + h * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
    y_rk10(T + k) = f(u_rk10(T + k) / theta);
  end
  
  for k = 1 : dt_pk10
    T = (t - 1) * dt_pk10;
    h = 1.0 / dt_pk10;
    u_p = u_pk10(T + k - 1) + h * du(u_pk10(T + k - 1), y_pk10(T + k - 1), i(t));
    u_pk10(T + k) = u_pk10(T + k - 1) + h * (du(u_pk10(T + k - 1), y_pk10(T + k - 1), i(t)) + du(u_p, f(u_p/theta), i(t))) / 2.0;
    y_pk10(T + k) = f(u_pk10(T + k) / theta);
  end
end

%% [sigmoidal model] differential models: plot results
figure();

subplot(5, 3, 1 : 2);
hold on; grid off; 
xlabel('t'); 
ylabel('i');
plot(1 : duration, i, 'Color', 'black');

subplot(5, 3, 4 : 5);
hold on; grid off; 
xlabel('t'); 
ylabel('y');
h_e10  = plot([1 : duration * dt_e10] ./ dt_e10,   y_e10,  ':k');
h_e100 = plot([1 : duration * dt_e100] ./ dt_e100, y_e100, '-.k');
h_pk10 = plot([1 : duration * dt_pk10] ./ dt_pk10, y_pk10, '--k');
h_rk10 = plot([1 : duration * dt_rk10] ./ dt_rk10, y_rk10, '-k');
y_rect = [103.0 0.3907 103.06 0.3910];
scatter(y_rect(1), y_rect(2), 'ok', 'fill');

subplot(5, 3, 6);
hold on; grid off; box on;
plot([1 : duration * dt_e10] ./ dt_e10,   y_e10,  ':k');
plot([1 : duration * dt_e100] ./ dt_e100, y_e100, '-.k');
plot([1 : duration * dt_pk10] ./ dt_pk10, y_pk10, '--k');
plot([1 : duration * dt_rk10] ./ dt_rk10, y_rk10, '-k');
xlim([y_rect(1) y_rect(3)]);
ylim([y_rect(2) y_rect(4)]);

subplot(5, 3, 7 : 8);
hold on; grid off; 
xlabel('t'); 
ylabel('\Delta y');
plot(1 : duration, abs(y_rk10([1 : duration] .* dt_rk10) - y_e10([1 : duration] .* dt_e10)),   ':k');
plot(1 : duration, abs(y_rk10([1 : duration] .* dt_rk10) - y_e100([1 : duration] .* dt_e100)), '-.k');
plot(1 : duration, abs(y_rk10([1 : duration] .* dt_rk10) - y_pk10([1 : duration] .* dt_pk10)), '--k');
dy_rect = [98 0 (110 - 98) (0.0009 - 0)];
rectangle('Position', dy_rect);

subplot(5, 3, 9);
hold on; grid off; box on;
plot(1 : duration, abs(y_rk10([1 : duration] .* dt_rk10) - y_e10([1 : duration] .* dt_e10)),   ':k');
plot(1 : duration, abs(y_rk10([1 : duration] .* dt_rk10) - y_e100([1 : duration] .* dt_e100)), '-.k');
plot(1 : duration, abs(y_rk10([1 : duration] .* dt_rk10) - y_pk10([1 : duration] .* dt_pk10)), '--k');
xlim([dy_rect(1) dy_rect(1) + dy_rect(3)]);
ylim([dy_rect(2) dy_rect(2) + dy_rect(4)]);

subplot(5, 3, 10 : 11);
hold on; grid off;
xlabel('t'); 
ylabel('u');
plot([1 : duration * dt_e10] ./ dt_e10,   u_e10,  ':k');
plot([1 : duration * dt_e100] ./ dt_e100, u_e100, '-.k');
plot([1 : duration * dt_pk10] ./ dt_pk10, u_pk10, '--k');
plot([1 : duration * dt_rk10] ./ dt_rk10, u_rk10, '-k');
y_rect = [103.0 2.554 103.06 2.559];
scatter(y_rect(1), y_rect(2), 'ok', 'fill');

subplot(5, 3, 12);
hold on; grid off; box on;
plot([1 : duration * dt_e10] ./ dt_e10,   u_e10,  ':k');
plot([1 : duration * dt_e100] ./ dt_e100, u_e100, '-.k');
plot([1 : duration * dt_pk10] ./ dt_pk10, u_pk10, '--k');
plot([1 : duration * dt_rk10] ./ dt_rk10, u_rk10, '-k');
xlim([y_rect(1) y_rect(3)]);
ylim([y_rect(2) y_rect(4)]);

subplot(5, 3, 13 : 14);
hold on; grid off; 
xlabel('t'); 
ylabel('\Delta u');
plot(1 : duration, abs(u_rk10([1 : duration] .* dt_rk10) - u_e10([1 : duration] .* dt_e10)),   ':k');
plot(1 : duration, abs(u_rk10([1 : duration] .* dt_rk10) - u_e100([1 : duration] .* dt_e100)), '-.k');
plot(1 : duration, abs(u_rk10([1 : duration] .* dt_rk10) - u_pk10([1 : duration] .* dt_pk10)), '--k');
du_rect = [98 0 (110 - 98) (0.003 - 0)];
rectangle('Position', du_rect);

subplot(5, 3, 15);
hold on; grid off; box on;
plot(1 : duration, abs(u_rk10([1 : duration] .* dt_rk10) - u_e10([1 : duration] .* dt_e10)),   ':k');
plot(1 : duration, abs(u_rk10([1 : duration] .* dt_rk10) - u_e100([1 : duration] .* dt_e100)), '-.k');
plot(1 : duration, abs(u_rk10([1 : duration] .* dt_rk10) - u_pk10([1 : duration] .* dt_pk10)), '--k');
xlim([du_rect(1) du_rect(1) + du_rect(3)]);
ylim([du_rect(2) du_rect(2) + du_rect(4)]);

legend([h_e10 h_e100 h_pk10 h_rk10], 'метод Эйлера (h=0.1)', 'метод Эйлера (h=0.01)', sprintf('%s\n%s', 'метод Эйлера с', 'пересчётом (h=0.1)'), 'метод Рунге-Кутты (h=0.1)');

%% [original model] differential models: model simulation
mu = 0.75;
threshold = 1.0;
alpha = 3.5;
theta = 1.0;

f = @(u) f_origin(u);
du = @(u, y, i) alpha .* y + i - threshold - mu .* u;

duration = 175;
i = zeros(1, duration);
i( 25 :  50) = 0.5;
i( 50 :  75) = 1.2;
i( 75 : 100) = 3.5;
i(100 : 125) = 0.5;
i(125 : 150) = 0.0;

dt_e10 = 10;
u_e10 = zeros(1, duration * dt_e10);
y_e10 = zeros(1, duration * dt_e10);

dt_e100 = 100;
u_e100 = zeros(1, duration * dt_e100);
y_e100 = zeros(1, duration * dt_e100);

dt_rk10 = 10;
u_rk10 = zeros(1, duration * dt_rk10);
y_rk10 = zeros(1, duration * dt_rk10);

dt_pk10 = 10;
u_pk10 = zeros(1, duration * dt_pk10);
y_pk10 = zeros(1, duration * dt_pk10);

for t = 2 : duration
  for k = 1 : dt_e10
    T = (t - 1) * dt_e10;
    u_e10(T + k) = u_e10(T + k - 1) + du(u_e10(T + k - 1), y_e10(T + k - 1), i(t)) / dt_e10;
    y_e10(T + k) = f(u_e10(T + k) / theta);
  end
  
  for k = 1 : dt_e100
    T = (t - 1) * dt_e100;
    u_e100(T + k) = u_e100(T + k - 1) + du(u_e100(T + k - 1), y_e100(T + k - 1), i(t)) / dt_e100;
    y_e100(T + k) = f(u_e100(T + k) / theta);
  end
  
  for k = 1 : dt_rk10
    T = (t - 1) * dt_rk10;
    h = 1.0 / dt_rk10;
    k1 = du(u_rk10(T + k - 1), y_rk10(T + k - 1), i(t));
    k2 = du(u_rk10(T + k - 1) + k1 * h / 2.0, f((u_rk10(T + k - 1) + k1 * h / 2.0) / theta), i(t));
    k3 = du(u_rk10(T + k - 1) + k2 * h / 2.0, f((u_rk10(T + k - 1) + k2 * h / 2.0) / theta), i(t));
    k4 = du(u_rk10(T + k - 1) + k3 * h, f((u_rk10(T + k - 1) + k3 * h) / theta), i(t));
    u_rk10(T + k) = u_rk10(T + k - 1) + h * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
    y_rk10(T + k) = f(u_rk10(T + k) / theta);
  end
  
  for k = 1 : dt_pk10
    T = (t - 1) * dt_pk10;
    h = 1.0 / dt_pk10;
    u_p = u_pk10(T + k - 1) + h * du(u_pk10(T + k - 1), y_pk10(T + k - 1), i(t));
    u_pk10(T + k) = u_pk10(T + k - 1) + h * (du(u_pk10(T + k - 1), y_pk10(T + k - 1), i(t)) + du(u_p, f(u_p/theta), i(t))) / 2.0;
    y_pk10(T + k) = f(u_pk10(T + k) / theta);
  end
end

%% [original model] differential models: plot results
figure();

subplot(5, 3, 1 : 2);
hold on; grid off; 
xlabel('t'); 
ylabel('i');
plot(1 : duration, i, 'Color', 'black');

subplot(5, 3, 4 : 5);
hold on; grid off; 
xlabel('t'); 
ylabel('y');
h_e10  = plot([1 : duration * dt_e10] ./ dt_e10,   y_e10,  ':k');
h_e100 = plot([1 : duration * dt_e100] ./ dt_e100, y_e100, '-.k');
h_pk10 = plot([1 : duration * dt_pk10] ./ dt_pk10, y_pk10, '--k');
h_rk10 = plot([1 : duration * dt_rk10] ./ dt_rk10, y_rk10, '-k');
y_rect = [102.93 0.68 103.04 0.6806];
scatter(y_rect(1), y_rect(2), 'ok', 'fill');

subplot(5, 3, 6);
hold on; grid off; box on;
plot([1 : duration * dt_e10] ./ dt_e10,   y_e10,  ':k');
plot([1 : duration * dt_e100] ./ dt_e100, y_e100, '-.k');
plot([1 : duration * dt_pk10] ./ dt_pk10, y_pk10, '--k');
plot([1 : duration * dt_rk10] ./ dt_rk10, y_rk10, '-k');
xlim([y_rect(1) y_rect(3)]);
ylim([y_rect(2) y_rect(4)]);

subplot(5, 3, 7 : 8);
hold on; grid off; 
xlabel('t'); 
ylabel('\Delta y');
plot(1 : duration, abs(y_rk10([1 : duration] .* dt_rk10) - y_e10([1 : duration] .* dt_e10)),   ':k');
plot(1 : duration, abs(y_rk10([1 : duration] .* dt_rk10) - y_e100([1 : duration] .* dt_e100)), '-.k');
plot(1 : duration, abs(y_rk10([1 : duration] .* dt_rk10) - y_pk10([1 : duration] .* dt_pk10)), '--k');
dy_rect = [98 0 (112 - 98) (0.002 - 0)];
rectangle('Position', dy_rect);

subplot(5, 3, 9);
hold on; grid off; box on;
plot(1 : duration, abs(y_rk10([1 : duration] .* dt_rk10) - y_e10([1 : duration] .* dt_e10)),   ':k');
plot(1 : duration, abs(y_rk10([1 : duration] .* dt_rk10) - y_e100([1 : duration] .* dt_e100)), '-.k');
plot(1 : duration, abs(y_rk10([1 : duration] .* dt_rk10) - y_pk10([1 : duration] .* dt_pk10)), '--k');
xlim([dy_rect(1) dy_rect(1) + dy_rect(3)]);
ylim([dy_rect(2) dy_rect(2) + dy_rect(4)]);

subplot(5, 3, 10 : 11);
hold on; grid off;
xlabel('t'); 
ylabel('u');
plot([1 : duration * dt_e10] ./ dt_e10,   u_e10,  ':k');
plot([1 : duration * dt_e100] ./ dt_e100, u_e100, '-.k');
plot([1 : duration * dt_pk10] ./ dt_pk10, u_pk10, '--k');
plot([1 : duration * dt_rk10] ./ dt_rk10, u_rk10, '-k');
y_rect = [103.96 2.371 104.005 2.375];
scatter(y_rect(1), y_rect(2), 'ok', 'fill');

subplot(5, 3, 12);
hold on; grid off; box on;
plot([1 : duration * dt_e10] ./ dt_e10,   u_e10,  ':k');
plot([1 : duration * dt_e100] ./ dt_e100, u_e100, '-.k');
plot([1 : duration * dt_pk10] ./ dt_pk10, u_pk10, '--k');
plot([1 : duration * dt_rk10] ./ dt_rk10, u_rk10, '-k');
xlim([y_rect(1) y_rect(3)]);
ylim([y_rect(2) y_rect(4)]);

subplot(5, 3, 13 : 14);
hold on; grid off; 
xlabel('t'); 
ylabel('\Delta u');
plot(1 : duration, abs(u_rk10([1 : duration] .* dt_rk10) - u_e10([1 : duration] .* dt_e10)),   ':k');
plot(1 : duration, abs(u_rk10([1 : duration] .* dt_rk10) - u_e100([1 : duration] .* dt_e100)), '-.k');
plot(1 : duration, abs(u_rk10([1 : duration] .* dt_rk10) - u_pk10([1 : duration] .* dt_pk10)), '--k');
du_rect = [98 0 (112 - 98) (0.009 - 0)];
rectangle('Position', du_rect);

subplot(5, 3, 15);
hold on; grid off; box on;
plot(1 : duration, abs(u_rk10([1 : duration] .* dt_rk10) - u_e10([1 : duration] .* dt_e10)),   ':k');
plot(1 : duration, abs(u_rk10([1 : duration] .* dt_rk10) - u_e100([1 : duration] .* dt_e100)), '-.k');
plot(1 : duration, abs(u_rk10([1 : duration] .* dt_rk10) - u_pk10([1 : duration] .* dt_pk10)), '--k');
xlim([du_rect(1) du_rect(1) + du_rect(3)]);
ylim([du_rect(2) du_rect(2) + du_rect(4)]);

legend([h_e10 h_e100 h_pk10 h_rk10], 'метод Эйлера (h=0.1)', 'метод Эйлера (h=0.01)', sprintf('%s\n%s', 'метод Эйлера с', 'пересчётом (h=0.1)'), 'метод Рунге-Кутты (h=0.1)');




