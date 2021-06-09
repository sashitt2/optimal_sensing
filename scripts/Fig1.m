% get the supercritical CGLE
nx = 220;
[xh,A,M,cd,U,mu0,mu2,cc,gam] = setup_supercritical(nx);
beta = 7;
q0 = 1/5/(sqrt(2*pi))*exp(-(xh-8).^2/2/5^2);
m = size(A,1);
Q = 0.01*eye(m);
R = 0.01;
rng(0);

tmax=200;  dt = 1;
tvec = 0:dt:tmax;
nt = length(tvec);
%%
%convert true system to discrete
dt = 1;
Ad = expm(A*dt);
Qd = Q*dt;
Rd = R*dt;
%%
q = zeros(size(q0,1),nt);
q(:,1) = q0;

for i=2:nt
        q(:,i) = Ad*q(:,i-1) + 0.01*randn(size(m,1));
end

% measurement noise
qm = q + 0.01*randn(size(q));

%%
% construct the sensor matrix and visualize its derivative w.r.t. xs
xs = 0.3128;

sigma = 0.25;
C1 = exp(-(xh-xs).^2/2/sigma^2)'*M;
C1 = C1 / max(C1);
dCdxs1 = C1.*(xh - xs)'/sigma^2;
sigma = 0.5;
C2 = exp(-(xh-xs).^2/2/sigma^2)'*M;
C2 = C2 / max(C2);
dCdxs2 = C2.*(xh - xs)'/sigma^2;
sigma = 1;
C3 = exp(-(xh-xs).^2/2/sigma^2)'*M;
C3 = C3 / max(C3);
dCdxs3 = C3.*(xh - xs)'/sigma^2;

figure('Position',[10 200 1000 800]);

plot(xh, C1', '-ob', xh, C2', '-^r', xh, C3, '-sk', ...
                'lineWidth', 3, 'MarkerSize', 15);

set(gca, 'FontSize', 30);
xlim([-10 10]);
ylabel('sensor profile', 'interpreter','latex');
xlabel('x', 'interpreter','latex')

leg = legend('$\Delta x/\sigma$ = 2.5', ...
             '$\Delta x/\sigma$ = 1.25', ...
             '$\Delta x/\sigma$ = 0.625');
set(leg,'Interpreter','latex','location','northeast');


ylim([-0.1, 1.1]);

h=gcf;                                                                                                                                     
set(h,'PaperOrientation','landscape');                                                                                                     
set(h,'PaperUnits','normalized');                                                                                                          
set(h,'PaperPosition', [0 0 1 1]);                                                                                                         
print(gcf, '-dpdf', 'figures/fig1_a.pdf')
%%
sigma = 0.25;

xs = -8:1:8;

J_inf1 = zeros(size(xs));
dJ_inf1 = zeros(size(xs));

istart = 5;
for i = istart:size(xs,2)
    
    C = exp(-(xh-xs(i)).^2/2/sigma^2)'*M;
    C = C / max(C);
    dCdxs = C.*(xh - xs(i))'/sigma^2;

    J_inf1(i) = getObjInfinite(Ad, qm(:,1), qm(:,2:end), Qd, Rd, C, eye(m));
    [~, ~, dJdC, ~] = ...
            getGradInfinite(Ad, qm(:,1), qm(:,2:end), Qd, Rd, C, eye(m));
        
    dJ_inf1(i) = real(sum(conj(dCdxs).*dJdC, 2));
    
    disp(['xs = ', num2str(xs(i)), ' has J = ', num2str(J_inf1(i)), ...
          ' and dJ = ', num2str(dJ_inf1(i))]);
end

sigma = 0.5;

xs = -8:1:8;

J_inf2 = zeros(size(xs));
dJ_inf2 = zeros(size(xs));

istart = 5;
for i = istart:size(xs,2)
    
    C = exp(-(xh-xs(i)).^2/2/sigma^2)'*M;
    C = C / max(C);
    dCdxs = C.*(xh - xs(i))'/sigma^2;

    J_inf2(i) = getObjInfinite(Ad, qm(:,1), qm(:,2:end), Qd, Rd, C, eye(m));
    [~, ~, dJdC, ~] = ...
            getGradInfinite(Ad, qm(:,1), qm(:,2:end), Qd, Rd, C, eye(m));
        
    dJ_inf2(i) = real(sum(conj(dCdxs).*dJdC, 2));
    
    disp(['xs = ', num2str(xs(i)), ' has J = ', num2str(J_inf2(i)), ...
          ' and dJ = ', num2str(dJ_inf2(i))]);
end


sigma = 1;

xs = -8:1:8;

J_inf3 = zeros(size(xs));
dJ_inf3 = zeros(size(xs));

istart = 5;
for i = istart:size(xs,2)
    
    C = exp(-(xh-xs(i)).^2/2/sigma^2)'*M;
    C = C / max(C);
    dCdxs = C.*(xh - xs(i))'/sigma^2;

    J_inf3(i) = getObjInfinite(Ad, qm(:,1), qm(:,2:end), Qd, Rd, C, eye(m));
    [~, ~, dJdC, ~] = ...
            getGradInfinite(Ad, qm(:,1), qm(:,2:end), Qd, Rd, C, eye(m));
        
    dJ_inf3(i) = real(sum(conj(dCdxs).*dJdC, 2));
    
    disp(['xs = ', num2str(xs(i)), ' has J = ', num2str(J_inf3(i)), ...
          ' and dJ = ', num2str(dJ_inf3(i))]);
end

%%
figure('Position',[10 200 1000 800]);

xs = -8:1:8;

istart = 5;
iend = size(xs,2);
dxs = 0.25;

const = 0.5 * norm(qm(:,1:end), 'fro')^2;

plot(xs(istart:iend), J_inf1(istart:iend)/const, ...
                'ok', ...
            'MarkerSize', 15, 'LineWidth', 3);
hold on

plot(xs(istart:iend), J_inf2(istart:iend)/const, ...
                '^k', ...
            'MarkerSize', 15, 'LineWidth', 3);

plot(xs(istart:iend), J_inf3(istart:iend)/const, ...
                'sk', ...
            'MarkerSize', 15, 'LineWidth', 3);

for i = istart:iend
    J = J_inf1(i)/const;
    dJ = dJ_inf1(i)/const;
    plot([xs(i) - dxs, xs(i) + dxs], [J - dJ*dxs, J + dJ*dxs], 'b', ...
            'LineWidth', 3);
end

for i = istart:iend
    J = J_inf2(i)/const;
    dJ = dJ_inf2(i)/const;
    plot([xs(i) - dxs, xs(i) + dxs], [J - dJ*dxs, J + dJ*dxs], 'r', ...
            'LineWidth', 3);
end


for i = istart:iend
    J = J_inf3(i)/const;
    dJ = dJ_inf3(i)/const;
    plot([xs(i) - dxs, xs(i) + dxs], [J - dJ*dxs, J + dJ*dxs], 'k', ...
            'LineWidth', 3);
end

hold off

ax = gca;
ax.XTick = [-5:2:10];

ylabel('normalized estimation error','interpreter','latex');
xlabel('sensor location','interpreter','latex')

leg = legend('$\Delta x/\sigma$ = 2.5', ...
             '$\Delta x/\sigma$ = 1.25', ...
             '$\Delta x/\sigma$ = 0.625');
set(leg,'Interpreter','latex','location','northeast');

set(gca, 'FontSize', 30);
h=gcf;                                                                                                                                     
set(h,'PaperOrientation','landscape');                                                                                                     
set(h,'PaperUnits','normalized');                                                                                                          
set(h,'PaperPosition', [0 0 1 1]);                                                                                                         

print(gcf, '-dpdf', 'figures/fig1_b.pdf')

%%
sigma = 0.25;
%xs = -1:0.2:1;
xs = 0:0.2:2;

resJ_inf1 = zeros(size(xs));
resdJ_inf1 = zeros(size(xs));

for i = 1:size(xs,2)
    
    C = exp(-(xh-xs(i)).^2/2/sigma^2)'*M;
    dCdxs = C.*(xh - xs(i))'/sigma^2;

    resJ_inf1(i) = getObjInfinite(Ad, qm(:,1), qm(:,2:end), Qd, Rd, C, eye(m));
    [~, ~, dJdC, ~] = ...
            getGradInfinite(Ad, qm(:,1), qm(:,2:end), Qd, Rd, C, eye(m));
        
    resdJ_inf1(i) = real(sum(conj(dCdxs).*dJdC, 2));
    
    disp(['xs = ', num2str(xs(i)), ' has J = ', num2str(resJ_inf1(i)), ...
          ' and dJ = ', num2str(resdJ_inf1(i))]);
end

sigma = 0.5;
xs = 0:0.2:2;


resJ_inf2 = zeros(size(xs));
resdJ_inf2 = zeros(size(xs));

for i = 1:size(xs,2)
    
    C = exp(-(xh-xs(i)).^2/2/sigma^2)'*M;
    dCdxs = C.*(xh - xs(i))'/sigma^2;

    resJ_inf2(i) = getObjInfinite(Ad, qm(:,1), qm(:,2:end), Qd, Rd, C, eye(m));
    [~, ~, dJdC, ~] = ...
            getGradInfinite(Ad, qm(:,1), qm(:,2:end), Qd, Rd, C, eye(m));
        
    resdJ_inf2(i) = real(sum(conj(dCdxs).*dJdC, 2));
    
    disp(['xs = ', num2str(xs(i)), ' has J = ', num2str(resJ_inf2(i)), ...
          ' and dJ = ', num2str(resdJ_inf2(i))]);
end


sigma = 1;
xs = 0:0.2:2;


resJ_inf3 = zeros(size(xs));
resdJ_inf3 = zeros(size(xs));

for i = 1:size(xs,2)
    
    C = exp(-(xh-xs(i)).^2/2/sigma^2)'*M;
    dCdxs = C.*(xh - xs(i))'/sigma^2;

    resJ_inf3(i) = getObjInfinite(Ad, qm(:,1), qm(:,2:end), Qd, Rd, C, eye(m));
    [~, ~, dJdC, ~] = ...
            getGradInfinite(Ad, qm(:,1), qm(:,2:end), Qd, Rd, C, eye(m));
        
    resdJ_inf3(i) = real(sum(conj(dCdxs).*dJdC, 2));
    
    disp(['xs = ', num2str(xs(i)), ' has J = ', num2str(resJ_inf3(i)), ...
          ' and dJ = ', num2str(resdJ_inf3(i))]);
end
%%
% resolved plot

figure('Position',[10 200 1000 800]);

xs = 0:0.2:2;

istart = 1;
iend = size(xs,2);
dxs = 0.1;

const = 0.5 * norm(qm(:,1:end), 'fro')^2;

plot(xs(istart:iend), resJ_inf1(istart:iend)/const, 'ok', ...
            'MarkerSize', 15, 'LineWidth', 3);

hold on

plot(xs(istart:iend), resJ_inf2(istart:iend)/const, '^k', ...
            'MarkerSize', 15, 'LineWidth', 3);

plot(xs(istart:iend), resJ_inf3(istart:iend)/const, 'sk', ...
            'MarkerSize', 15, 'LineWidth', 3);        

for i = istart:iend
    J = resJ_inf1(i)/const;
    dJ = resdJ_inf1(i)/const;
    plot([xs(i) - dxs, xs(i) + dxs], [J - dJ*dxs, J + dJ*dxs], 'b', ...
            'LineWidth', 3);
end

for i = istart:iend
    J = resJ_inf2(i)/const;
    dJ = resdJ_inf2(i)/const;
    plot([xs(i) - dxs, xs(i) + dxs], [J - dJ*dxs, J + dJ*dxs], 'r', ...
            'LineWidth', 3);
end

for i = istart:iend
    J = resJ_inf3(i)/const;
    dJ = resdJ_inf3(i)/const;
    plot([xs(i) - dxs, xs(i) + dxs], [J - dJ*dxs, J + dJ*dxs], 'k', ...
            'LineWidth', 3);
end


hold off

ylabel('normalized estimation error','interpreter','latex');
xlabel('sensor location','interpreter','latex')

leg = legend('$\Delta x/\sigma$ = 2.5', ...
             '$\Delta x/\sigma$ = 1.25', ...
             '$\Delta x/\sigma$ = 0.625');
set(leg,'Interpreter','latex','location','northeast');

xlim([-0.2 2.2])

set(gca, 'FontSize', 30);
h=gcf;                                                                                                                                     
set(h,'PaperOrientation','landscape');                                                                                                     
set(h,'PaperUnits','normalized');                                                                                                          
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', 'figures/fig1_c.pdf')