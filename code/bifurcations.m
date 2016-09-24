%% common part for all blocks

k12 = 2.6 * logsig(38.73 / 3.5) - 0.49;
u01 = 2.6 * logsig(-23.22);
u12 = 2.6 * logsig(38.73 / 3.5);
f_origin = @(u) heaviside_restricted(u - u12) .* (1.0 - 0.35 ./ (u - k12)) + ...
                (heaviside_restricted(u - u01) - heaviside_restricted(u - u12)) .* (0.1935 + log(u ./ (2.6 - u)) ./ 120.0);
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


%% [original model] bifurcation: equlibrium intervals
figure();
hold on;
grid on;

x_bound = 1.0 / 3.5;
x = 0.01 : 0.01 : x_bound;
k = 1.0 ./ (1.0 + exp((1.0 ./ x - 5.2) ./ 0.23));
y = 2.6 .* k .* (1.0 - k) ./ (0.23 .* x.^2);

plot(x, y);

x = x_bound : 0.01 : 0.95;
y = 0.35 ./ (x - 1.0).^2;

plot(x, y);

mu = 0.75;
threshold = 1.0;
theta = 1.0;
alpha = 6.5;

plot([0.0 1.0], [alpha./(mu.*theta) alpha./(mu.*theta)]);


%% [original model] bifurcation: i vs alpha diagram
figure();
hold on;
grid on;
xlabel('i');
ylabel('\alpha');

mu = 0.75;
threshold = 1.0;
theta = 1.0;

%alpha = 1.4 * mu * theta / 3.0 : 0.1 : 25;
alpha = 1.4 * mu * theta / 3.0  .* [1 2 3];
f1 = (1.0 + sqrt(-1.4 .* mu .* theta ./ alpha + 3.0)) ./ 2.0;
i1 = alpha .* f1 - threshold - mu .* theta .* (2.11 + 0.35 ./ (1.0 - f1));
f2 = (1.0 - sqrt(-1.4 .* mu .* theta ./ alpha + 3.0)) ./ 2.0;
i2 = alpha .* f2 - threshold - mu .* theta .* (2.11 + 0.35 ./ (1.0 - f2));

plot(i1, alpha, 'Color', 'green');
plot(i2, alpha, 'Color', 'blue');










