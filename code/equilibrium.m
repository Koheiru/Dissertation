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

%% [sigmoidal model] equilibrium points for soft and hard states
fh = figure();
figure_adjust(fh, [17.5 10.0]);

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

figure_subplot(2, 2, 1);
hold on; grid off; box on;
plot(I_soft, y, '-k');
plot([i_soft i_soft], [0 1], '--k');
xlabel('i');
ylabel('y');
xlim([-1 5]);
ylim([-0.02 1.02]);

figure_subplot(2, 2, 3);
hold on; grid off; box on;
plot(y, F_soft, '-k');
plot([0 1], [0 0], '--k');
xlabel('y');
ylabel('F(y)');
ylim([-3.0 3.0])

figure_subplot(2, 2, 2);
hold on; grid off; box on;
plot(I_hard, y, '-k');
plot([i_hard i_hard], [0 1], '--k');
xlabel('i');
ylabel('y');
xlim([-1 4]);
ylim([-0.02 1.02]);

figure_subplot(2, 2, 4);
hold on; grid off; box on;
plot(y, F_hard, '-k');
plot([0 1], [0 0], '--k');
xlabel('y');
ylabel('F(y)');
ylim([-3.0 3.0])


%% [sigmoidal model] equilibrium conditions
fh = figure();
figure_adjust(fh, [17.5 12]);

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

figure_subplot(2, 3, 1);
hold on; grid off; box on;
plot(y, F_soft, '-k');
plot([0 1], [0 0], '--k');
xlabel('y');
ylabel('F');
ylim([-3 3]);

figure_subplot(2, 3, 2);
hold on; grid off; box on;
plot(y, F_zero, '-k');
plot([0 1], [0 0], '--k');
xlabel('y');
ylabel('F');
ylim([-3 3]);

figure_subplot(2, 3, 3);
hold on; grid off; box on;
plot(y, F_hard, '-k');
plot([0 1], [0 0], '--k');
xlabel('y');
ylabel('F');
ylim([-3 3]);

y = 0.1 : 0.01 : 0.9;
dF_zero = alpha_zero - mu .* theta .* (1 ./ y + 1.0 ./ (1.0 - y));
dF_soft = alpha_soft - mu .* theta .* (1 ./ y + 1.0 ./ (1.0 - y));
dF_hard = alpha_hard - mu .* theta .* (1 ./ y + 1.0 ./ (1.0 - y));

figure_subplot(2, 3, 4);
hold on; grid off; box on;
plot(y, dF_soft, '-k');
plot([0 1], [0 0], '--k');
xlabel('y');
ylabel('F''');
ylim([-3 3]);

figure_subplot(2, 3, 5);
hold on; grid off; box on;
plot(y, dF_zero, '-k');
plot([0 1], [0 0], '--k');
xlabel('y');
ylabel('F''');
ylim([-3 3]);

figure_subplot(2, 3, 6);
hold on; grid off; box on;
plot(y, dF_hard, '-k');
plot([0 1], [0 0], '--k');
xlabel('y');
ylabel('F''');
ylim([-3 3]);


%% [sigmoidal model] equilibrium bound points
fh = figure();
figure_adjust(fh, [17.5 6.5]);

mu = 0.75;
threshold = 1.0;
theta     = 1.0;
i         = 1.0;
alpha_zero = 4.0 * mu * theta;
alpha_hard = alpha_zero * 1.5;

f_zero = 0.5;
f_hard = (1.0 - sqrt(1.0 - 4.0 .* mu .* theta ./ alpha_hard)) ./ 2.0;
i_zero = -(alpha_zero .* f_zero - threshold - mu .* theta .* s_sigm(f_zero));
i_hard = -(alpha_hard .* f_hard - threshold - mu .* theta .* s_sigm(f_hard));

y = 0.001 : 0.001 : 0.999;
F_zero = i_zero + alpha_zero .* y - threshold - mu .* theta .* s_sigm(y);
F_hard = i_hard + alpha_hard .* y - threshold - mu .* theta .* s_sigm(y);

figure_subplot(1, 2, 1);
hold on; grid off; box on;
plot(y, F_zero, '-k');
plot([0 1], [0 0], '--k');
xlabel('y');
ylabel('F');
xlim([0.46 0.54]);
ylim([-0.0002 0.0002]);
%set(gca, 'YTickMode','manual');
%set(gca, 'YTickLabel',num2str(get(gca,'YTick')'));

figure_subplot(1,2 , 2);
hold on; grid off; box on;
plot(y, F_hard, '-k');
plot([0 1], [0 0], '--k');
xlabel('y');
ylabel('F');
xlim([0.1 0.4]);
ylim([-0.1 0.3]);
%set(gca, 'YTickMode','manual');
%set(gca, 'YTickLabel',num2str(get(gca,'YTick')'));


%% [original model] equilibrium points for soft and hard states + special uninteresting state
fh = figure();
figure_adjust(fh, [17.5 20.5]);

mu = 0.75;
threshold = 1.0;
theta = 1.0;
alpha = [0.3 0.5  2.5 5.5  60.0];
i     = [0.0 0.0 2.0 0.5 -20.0];

y = [0.00 : 0.01 : 0.99, 0.991 : 0.001 : 0.999];
F = bsxfun(@plus, i', alpha' * y - threshold - repmat(mu .* theta .* s_origin(y), length(alpha), 1));
I = -(alpha' * y - threshold - repmat(mu .* theta .* s_origin(y), length(alpha), 1));

%----------------------
figure_subplot(5, 2, 1);
hold on; grid off; box on;
plot(I(1,:), y(:), '-k'); 
plot([-100 threshold], [0 0], '-k');
plot([i(1) i(1)], [0 1], '--k');
xlabel('i');
ylabel('y');
xlim([-0.5 3.5]);
ylim([-0.02 1.02]);

figure_subplot(5, 2, 2);
hold on; grid off; box on;
plot(y(:), F(1,:), '-k');
plot([0 0 ], [F(1,1) 100], '-k');
plot([0 1], [0 0], '--k');
xlabel('y');
ylabel('F(y)');
xlim([-0.02 1.02]);
ylim([-3.0 3.0])

%----------------------
rect_I = [0.6 -0.02; 1.4 0.22];
rect_F = [-0.02 -0.25; 0.22 0.25];

figure_subplot(5, 2, 3);
hold on; grid off; box on;
plot(I(2,:), y(:), '-k');
plot([-100 threshold], [0 0], '-k');
plot([i(2) i(2)], [0 1], '--k');
rectangle('Position', [rect_I(1,:), rect_I(2,:) - rect_I(1,:)]);
xlabel('i');
ylabel('y');
xlim([-0.5 3.5]);
ylim([-0.02 1.02]);

figure_subplot(5, 2, 4);
hold on; grid off; box on;
plot(y(:), F(2,:), '-k');
plot([0 0 ], [F(2,1) 100], '-k');
plot([0 1], [0 0], '--k');
rectangle('Position', [rect_F(1,:), rect_F(2,:) - rect_F(1,:)]);
xlabel('y');
ylabel('F(y)');
xlim([-0.02 1.02]);
ylim([-3.0 3.0])

%----------------------
figure_subplot(5, 2, 5);
hold on; grid off; box on;
plot(I(3,:), y(:), '-k');
plot([-100 threshold], [0 0], '-k');
plot([i(3) i(3)], [0 1], '--k');
xlabel('i');
ylabel('y');
xlim([-0.5 3.5]);
ylim([-0.02 1.02]);

figure_subplot(5, 2, 6);
hold on; grid off; box on;
plot(y(:), F(3,:), '-k');
plot([0 0 ], [F(3,1) 100], '-k');
plot([0 1], [0 0], '--k');
xlabel('y');
ylabel('F(y)');
xlim([-0.02 1.02]);
ylim([-2.0 2.0])

%----------------------
figure_subplot(5, 2, 7);
hold on; grid off; box on;
plot(I(4,:), y(:), '-k');
plot([-100 threshold], [0 0], '-k');
plot([i(4) i(4)], [0 1], '--k');
xlabel('i');
ylabel('y');
xlim([-0.5 3.5]);
ylim([-0.02 1.02]);

figure_subplot(5, 2, 8);
hold on; grid off; box on;
plot(y(:), F(4,:), '-k');
plot([0 0 ], [F(4,1) 100], '-k');
plot([0 1], [0 0], '--k');
xlabel('y');
ylabel('F(y)');
xlim([-0.02 1.02]);
ylim([-2.0 2.0])

%----------------------
figure_subplot(5, 2, 9);
hold on; grid off; box on;
plot(I(5,:), y(:), '-k');
plot([-100 threshold], [0 0], '-k');
plot([i(5) i(5)], [0 1], '--k');
xlabel('i');
ylabel('y');
xlim([-50   3.5]);
ylim([-0.02 1.02]);

figure_subplot(5, 2, 10);
hold on; grid off; box on;
plot(y(:), F(5,:), '-k');
plot([0 0 ], [F(5,1) 100], '-k');
plot([0 1], [0 0], '--k');
xlabel('y');
ylabel('F(y)');
xlim([-0.02 1.02]);
ylim([-30.0 30.0])

%% scaling some subplots
fh = figure();
figure_adjust(fh, [5.5 3.0]);

plot(I(2,:), y(:), '-k'); hold on; 
plot([-100 threshold], [0 0], '-k'); hold on; 
plot([i(2) i(2)], [0 1], '--k'); hold on; 
grid off; box on;
xlabel('i');
ylabel('y');
xlim(rect_I(:,1)');
ylim(rect_I(:,2)');

%----------------------
fh = figure();
figure_adjust(fh, [5.5 3.0]);

plot(y(:), F(2,:), '-k'); hold on;
plot([0 0], [F(2,1), 100], '-k'); hold on;
plot([0 1], [0 0], '--k'); hold on;
grid off; box on;
xlabel('y');
ylabel('F(y)');
xlim(rect_F(:,1)');
ylim(rect_F(:,2)');


%% [original model] equilibrium conditions
fh = figure();
figure_adjust(fh, [17.5 12]);

mu = 0.75;
threshold = 1.0;
theta     = 1.0;
alpha     = 75.0;

y = 0.01 : 0.001 : y12;
x = exp(-120 .* (y - 0.1935));
dF = alpha - mu .* theta .* 312 .* x ./ ((1.0 + x) .^ 2);

plot(y, dF); hold on;
plot([0 0], [min(dF) max(dF)]);
plot([y12 y12], [min(dF) max(dF)]);
grid off; box on;
xlabel('y');
ylabel('F^{''}(y)');

% mu = 0.75;
% threshold = 1.0;
% theta     = 1.0;
% i         = 1.0;
% alpha_zero = 4.0 * mu * theta;
% alpha_soft = alpha_zero / 1.5;
% alpha_hard = alpha_zero * 1.5;
%
% y = 0.001 : 0.001 : 0.999;
% F_zero = i + alpha_zero .* y - threshold - mu .* theta .* s_sigm(y);
% F_soft = i + alpha_soft .* y - threshold - mu .* theta .* s_sigm(y);
% F_hard = i + alpha_hard .* y - threshold - mu .* theta .* s_sigm(y);
% 
% figure_subplot(2, 3, 1);
% hold on; grid off; box on;
% plot(y, F_soft, '-k');
% plot([0 1], [0 0], '--k');
% xlabel('y');
% ylabel('F');
% ylim([-3 3]);
% 
% figure_subplot(2, 3, 2);
% hold on; grid off; box on;
% plot(y, F_zero, '-k');
% plot([0 1], [0 0], '--k');
% xlabel('y');
% ylabel('F');
% ylim([-3 3]);
% 
% figure_subplot(2, 3, 3);
% hold on; grid off; box on;
% plot(y, F_hard, '-k');
% plot([0 1], [0 0], '--k');
% xlabel('y');
% ylabel('F');
% ylim([-3 3]);
% 
% y = 0.1 : 0.01 : 0.9;
% dF_zero = alpha_zero - mu .* theta .* (1 ./ y + 1.0 ./ (1.0 - y));
% dF_soft = alpha_soft - mu .* theta .* (1 ./ y + 1.0 ./ (1.0 - y));
% dF_hard = alpha_hard - mu .* theta .* (1 ./ y + 1.0 ./ (1.0 - y));
% 
% figure_subplot(2, 3, 4);
% hold on; grid off; box on;
% plot(y, dF_soft, '-k');
% plot([0 1], [0 0], '--k');
% xlabel('y');
% ylabel('F''');
% ylim([-3 3]);
% 
% figure_subplot(2, 3, 5);
% hold on; grid off; box on;
% plot(y, dF_zero, '-k');
% plot([0 1], [0 0], '--k');
% xlabel('y');
% ylabel('F''');
% ylim([-3 3]);
% 
% figure_subplot(2, 3, 6);
% hold on; grid off; box on;
% plot(y, dF_hard, '-k');
% plot([0 1], [0 0], '--k');
% xlabel('y');
% ylabel('F''');
% ylim([-3 3]);




