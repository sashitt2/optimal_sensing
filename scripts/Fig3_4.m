% get the supercritical CGLE
nx = 220;
[xh,A,M,cd,U,mu0,mu2,cc,gam] = setup_supercritical(nx);
beta = 7;
q0 = 1/5/(sqrt(2*pi))*exp(-(xh-8).^2/2/5^2);
m = size(A,1);
q_noise = 0.01;
r_noise = 0.01;
rng(0);

tmax=200;  dt = 1;
tvec = 0:dt:tmax;
nt = length(tvec);

trueColor = [0.8500 0.3250 0.0980];
dmdColor = [0.4940 0.1840 0.5560];
%%
%convert true system to discrete
dt = 1;
Ad = expm(A*dt);
%%
q = zeros(size(q0,1),nt);
q(:,1) = q0;

for i=2:nt
        q(:,i) = Ad*q(:,i-1) + 0.01*randn(size(m,1));
end

% measurement noise
qm = q + 0.01*randn(size(q));

% plot eigenvalues of the discrete operator
eAd = eig(Ad);
%%
% prediction of the flow
nt_long = 400;
qp = zeros(size(q0,1),nt_long);
qp(:,1) = q0;

for i=2:nt_long
        qp(:,i) = Ad*qp(:,i-1) + 0.01*randn(size(m,1));
end

% measurement noise
qpm = qp + 0.01*randn(size(qp));
%%
% data has been generated. Now we look at the optimal sensor location for
% varying number of sensor locations

nsList = [1, 2, 3, 4, 5];
%nsList = [3, 4, 5];
%nsList = 1;

optimalSensors = cell(numel(nsList), 2);
timeTaken = zeros(numel(nsList), 2);
objConv = cell(numel(nsList),2);
reconError = cell(numel(nsList),2);
predError = cell(numel(nsList),2);

Qd = q_noise*(eye(m));

r = 10;
[Ldmd, Mdmd, ~] = get_low_rank(qm(:,2:end),qm(:,1:end-1),r,'dmd');

update = [1, 0 , 0];

for iter = 1:numel(nsList)

    ns = nsList(iter);
    
    %xs0 = zeros(ns,1);
    xs0 = randn(ns, 1);
    
    Qdmd = Ldmd'*Qd*Ldmd;
    Rd = r_noise*(eye(ns));
    
    disp(['true dynamics optimal sensor location -> ', num2str(iter)]);
    tic
    [xsol, ~, ~, objList] = gradDescent(xs0, Ad, eye(m), xh, qm, Qd, Rd, ...
        update, 50, 2);
    timeTrue = toc;
    
    disp(['DMD dynamics optimal sensor location -> ', num2str(iter)]);
    tic
    [xsolDMD, ~, ~, objListDMD] = gradDescent(xs0, Mdmd, Ldmd, ...
                                        xh, qm, Qdmd, Rd, update, 50, 2);
    timeDMD = toc;
    
    timeTaken(iter,1) = timeTrue;
    timeTaken(iter,2) = timeDMD;
    
    optimalSensors{iter,1} = xsol;
    optimalSensors{iter,2} = xsolDMD;

    objConv{iter,1} = objList;
    objConv{iter,2} = objListDMD;
    
    reconError{iter,1} = evalSensorLocation(xsol,Ad,eye(m),...
                                            xh, qm, Qd, Rd);
    reconError{iter,2} = evalSensorLocation(xsolDMD,Mdmd, Ldmd,...
                                        xh, qm, Qdmd, Rd);
                                    
    predError{iter,1} = evalSensorLocation(xsol,Ad,eye(m),...
                                            xh, qpm, Qd, Rd);
    predError{iter,2} = evalSensorLocation(xsolDMD,Mdmd, Ldmd,...
                                        xh, qpm, Qdmd, Rd);                                    
end
%%
% sensor position plot

figure('Position',[10 200 1000 800]);
ywidth = 1;
yvalue = 1;
hold on
ymax = (numel(nsList))*ywidth;
xAmp = [-8.6 8.6 8.6 -8.6];
yAmp = [1 1 ymax+1 ymax+1];
patch( xAmp, yAmp, 'k','FaceAlpha', .1, 'EdgeColor', 'none');
for i = 1:numel(nsList)
    trueSensors = optimalSensors{i,1};
    
    for j = 1:numel(trueSensors)
        plot(trueSensors(j), yvalue, 'o', 'MarkerSize', 20, 'LineWidth', 3, ...
                        'MarkerEdgeColor', trueColor);
    end
    
    dmdSensors = optimalSensors{i,2};
    for j = 1:numel(trueSensors)
        plot(dmdSensors(j), yvalue, 'x', 'MarkerSize', 20, 'LineWidth', 3, ...
                        'MarkerEdgeColor', dmdColor);
    end
    
    yvalue = yvalue + ywidth;
end
hold off
set(gca, 'FontSize', 30)
%xlim([-12, 12]);
ylim([1, ymax+1]);
xlabel('sensor locations','interpreter','latex');
set(gca,'ytick',nsList)
ylabel('number of sensors','interpreter','latex');
line_vars = findall(gca);
leg = legend(line_vars([31,30]), {'full-order', 'DMD'});
set(leg,'Interpreter','latex','location','northeast');
box on

h = gcf;
set(h,'PaperOrientation','landscape');                                                                                                     
set(h,'PaperUnits','normalized');                                                                                                          
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', 'figures/fig3_a.pdf')


%%
figure('Position',[10 200 1000 800]);

reconErrorMat = cell2mat(reconError);

const = 0.5 * norm(qm, 'fro')^2;

plot(nsList, reconErrorMat(:,1) / const, '-o', 'MarkerSize', 15, 'LineWidth', 3, ...
                'Color', trueColor)
hold on
plot(nsList, reconErrorMat(:,2) / const ,'-x', 'MarkerSize', 15, 'LineWidth', 3, ...
                'Color', dmdColor)
set(gca, 'FontSize', 30)
leg = legend('full-order', 'DMD');

set(leg,'Interpreter','latex','location','northeast');
xlabel('number of sensors','interpreter','latex');
ylabel('normalized estimation error','interpreter','latex');
%%
nsList = [1, 2, 3, 4, 5];
figure('Position',[10 200 1000 800]);

semilogy(nsList, timeTaken(:,1), '-ok', 'MarkerSize', 15, 'LineWidth', 3, ...
                'Color', trueColor)
hold on
semilogy(nsList, timeTaken(:,2) ,'-xr', 'MarkerSize', 15, 'LineWidth', 3, ...
                'Color', dmdColor)
set(gca, 'FontSize', 30)
leg = legend('full-order', 'DMD');
set(leg,'Interpreter','latex','location','northeast');
xlabel('number of sensors','interpreter','latex');
ylabel('time taken (s)','interpreter','latex');

h = gcf;
set(h,'PaperOrientation','landscape');                                                                                                     
set(h,'PaperUnits','normalized');                                                                                                          
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', 'figures/fig3_b.pdf')
%%
% prediction and reconstruction analysis

nsList = [1, 2, 3, 4, 5];
%ntlong_list = [50, 100, 150, 200, 250, 300, 350, 400];
ntlong_list = [100, 200, 300, 400];

niter = 5;

trueError = zeros(numel(nsList), numel(ntlong_list), niter);
dmdError = zeros(numel(nsList), numel(ntlong_list), niter);

for ns_iter = 1:numel(nsList)
    
    ns = nsList(ns_iter);
    
    % Q and R matrices
    Qdmd = Ldmd'*Qd*Ldmd;
    Rd = r_noise*(eye(ns));
        
    xsol = optimalSensors{ns_iter,1};
    xsolDMD = optimalSensors{ns_iter,2};
    
    for idx = 1:numel(ntlong_list)

        ntlong = ntlong_list(idx);

        for iter = 1:niter

            disp([num2str(ns),'---',num2str(ntlong),'---',num2str(iter)]);
            
            % simulate noisy data
            qp = zeros(size(q0,1),ntlong);
            qp(:,1) = q0;

            for i=2:ntlong
                    qp(:,i) = Ad*qp(:,i-1) + q_noise*randn(size(m,1));
            end

            % measurement noise
            qpm = qp + r_noise*randn(size(qp));
                        
%             const = 0.5 * norm(qpm, 'fro')^2;
            % find prediction error
            trueError(ns_iter, idx, iter) = ...
                                evalSensorLocation(xsol,Ad,eye(m),...
                                                xh, qpm, Qd, Rd);
            dmdError(ns_iter, idx, iter) = ...
                                evalSensorLocation(xsolDMD,Mdmd, Ldmd,...
                                            xh, qpm, Qdmd, Rd);
        end
    end
end
%%

disp('---------------')
for iter = 1:niter
    disp([num2str(trueError(end,1,iter)),'---',num2str(dmdError(end,1,iter))])
end

%%
% prediction and reconstruction analysis
rng(0);
nsList = [1, 2, 3, 4, 5];

%ntlong_list = [50, 100, 150, 200, 250, 300, 350, 400];
ntlong_list = [100, 200, 300, 400];

niter = 20;

trueError = zeros(numel(nsList), numel(ntlong_list), niter);
dmdError = zeros(numel(nsList), numel(ntlong_list), niter);

q_noise = 0.01;
r_noise = 0.01;
Qd = q_noise*(eye(m));

for ns_iter = 1:numel(nsList)
    
    ns = nsList(ns_iter);
    
    % Q and R matrices
    Qdmd = Ldmd'*Qd*Ldmd;
    Rd = r_noise*(eye(ns));
        
    xsol = optimalSensors{ns_iter,1};
    xsolDMD = optimalSensors{ns_iter,2};
    
    for idx = 1:numel(ntlong_list)

        ntlong = ntlong_list(idx);

        for iter = 1:niter

            disp([num2str(ns),'---',num2str(ntlong),'---',num2str(iter)]);
            
            % simulate noisy data
            q0 = 1/5/(sqrt(2*pi))*exp(-(xh-8 + randn).^2/2/5^2);

            qp = zeros(size(q0,1),ntlong);
            qp(:,1) = q0;

            for i=2:ntlong
                    qp(:,i) = Ad*qp(:,i-1) + q_noise*randn(size(m,1));
            end

            % measurement noise
            qpm = qp + r_noise*randn(size(qp));
                        
%             const = 0.5 * norm(qpm, 'fro')^2;
            % find prediction error
            trueError(ns_iter, idx, iter) = ...
                                evalSensorLocation(xsol,Ad,eye(m),...
                                                xh, qpm, Qd, Rd);
            dmdError(ns_iter, idx, iter) = ...
                                evalSensorLocation(xsolDMD,Mdmd, Ldmd,...
                                            xh, qpm, Qdmd, Rd);
        end
    end
end
%%

trueColor = [0.8500 0.3250 0.0980];
dmdColor = [0.4940 0.1840 0.5560];

for ns_iter = 1:1%numel(nsList)
    figure('Position',[10 200 1000 800]);
    
    data1 = zeros(niter, numel(ntlong_list));
    for i = 1:numel(ntlong_list)
        data1(:,i) = trueError(ns_iter,i,:);
    end
    
    data2 = zeros(niter, numel(ntlong_list));
    for i = 1:numel(ntlong_list)
        data2(:,i) = dmdError(ns_iter,i,:);
    end    

    b = boxplot(data1./data2, ...
            'Labels', {'100', '200', '300', '400'}, ...
            'BoxStyle', 'outline', ...
            'Colors', 'k', 'Symbol','w+');
    set(b,{'linew'},{2})
    
%     hold on
% 
%     b = boxplot(data2, ...
%             'Labels', {'100', '200', '300', '400'}, ...
%             'BoxStyle', 'outline', ...
%             'Colors', dmdColor);
% 	set(b,{'linew'},{2})    

    hold on

    jitter = 0.05;
    
    x1 = ones(length(data1),1);
    scatter(x1 + jitter*randn(size(data1(:,1))),data1(:,1)./data2(:,1), 20, ...
                'filled', 'k', 'o', 'MarkerEdgeColor', 'k');
    scatter(x1*2 + jitter*randn(size(data1(:,1))),data1(:,2)./data2(:,2), 20, 'filled', 'k', 'o', 'MarkerEdgeColor', 'k');
    scatter(x1*3 + jitter*randn(size(data1(:,1))),data1(:,3)./data2(:,3), 20, 'filled', 'k', 'o', 'MarkerEdgeColor', 'k');
    scatter(x1*4 + jitter*randn(size(data1(:,1))),data1(:,4)./data2(:,4), 20, 'filled', 'k', 'o', 'MarkerEdgeColor', 'k');

    hold on
    
    plot(0:5, ones(6,1), '-.k', 'linewidth', 1);
    
    hold off

    
%     box_vars = findall(gca,'Tag','Box');
%     leg = legend(box_vars([5,1]), {'full-order', 'DMD'});
%     set(leg,'Interpreter','latex','location','northwest');
    
    set(gca, 'FontSize', 30)
    %legend('true', 'dmd')
    xlabel('time horizon','interpreter','latex');
    ylabel('relative estimation error','interpreter','latex');
    
end

h = gcf;
set(h,'PaperOrientation','landscape');                                                                                                     
set(h,'PaperUnits','normalized');                                                                                                          
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', 'figures/fig4_a.pdf')

%%
trueColor = [0.8500 0.3250 0.0980];
dmdColor = [0.4940 0.1840 0.5560];

for ns_iter = numel(nsList):numel(nsList)
    figure('Position',[10 200 1000 800]);
    
    data1 = zeros(niter, numel(ntlong_list));
    for i = 1:numel(ntlong_list)
        data1(:,i) = trueError(ns_iter,i,:);
    end
    
    data2 = zeros(niter, numel(ntlong_list));
    for i = 1:numel(ntlong_list)
        data2(:,i) = dmdError(ns_iter,i,:);
    end    

    b = boxplot(data1./data2, ...
            'Labels', {'100', '200', '300', '400'}, ...
            'BoxStyle', 'outline', ...
            'Colors', 'k', 'Symbol','w+');
    set(b,{'linew'},{2})
    
    hold on

    jitter = 0.05;
    
    x1 = ones(length(data1),1);
    scatter(x1 + jitter*randn(size(data1(:,1))),data1(:,1)./data2(:,1), 20, ...
                'filled', 'k', 'o', 'MarkerEdgeColor', 'k');
    scatter(x1*2 + jitter*randn(size(data1(:,1))),data1(:,2)./data2(:,2), 20, 'filled', 'k', 'o', 'MarkerEdgeColor', 'k');
    scatter(x1*3 + jitter*randn(size(data1(:,1))),data1(:,3)./data2(:,3), 20, 'filled', 'k', 'o', 'MarkerEdgeColor', 'k');
    scatter(x1*4 + jitter*randn(size(data1(:,1))),data1(:,4)./data2(:,4), 20, 'filled', 'k', 'o', 'MarkerEdgeColor', 'k');

    hold on
    
    plot(0:5, ones(6,1), '-.k', 'linewidth', 1);
    
    hold off

    
%     box_vars = findall(gca,'Tag','Box');
%     leg = legend(box_vars([5,1]), {'full-order', 'DMD'});
%     set(leg,'Interpreter','latex','location','northwest');
    
    set(gca, 'FontSize', 30)
    %legend('true', 'dmd')
    xlabel('time horizon','interpreter','latex');
    ylabel('relative estimation error','interpreter','latex');
    
end

h = gcf;
set(h,'PaperOrientation','landscape');                                                                                                     
set(h,'PaperUnits','normalized');                                                                                                          
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', 'figures/fig4_b.pdf')

%%
% noise level analysis
% horizon t = 200
nt = 201;
% prediction and reconstruction analysis
rng(0);
%nsList = [1, 2, 3, 4, 5];
nsList = [1];

%ntlong_list = [50, 100, 150, 200, 250, 300, 350, 400];
%ntlong_list = [100, 200, 300, 400];

noise_list = [0.005, 0.01, 0.02, 0.05];

niter = 20;

trueError_noise = zeros(numel(nsList), numel(noise_list), niter);
dmdError_noise = zeros(numel(nsList), numel(noise_list), niter);

for ns_iter = 1:numel(nsList)
    
    ns = nsList(ns_iter);
            
    xsol_noise = optimalSensors{ns_iter,1};
    xsolDMD_noise = optimalSensors{ns_iter,2};
    
    for idx = 1:numel(noise_list)

        curr_q_noise = noise_list(idx);
        curr_r_noise = noise_list(idx);
        Qd = curr_q_noise*(eye(m));

        % Q and R matrices
        Qdmd = Ldmd'*Qd*Ldmd;
        Rd = curr_r_noise*(eye(ns));        
        
        for iter = 1:niter

            disp([num2str(ns),'---',num2str(curr_q_noise),'---',num2str(iter)]);
            
            % simulate noisy data
            q0 = 1/5/(sqrt(2*pi))*exp(-(xh-8 + randn).^2/2/5^2);

            qp = zeros(size(q0,1),nt);
            qp(:,1) = q0;

            for i=2:nt
                    qp(:,i) = Ad*qp(:,i-1) + curr_q_noise*randn(size(m,1));
            end

            % measurement noise
            qpm = qp + curr_r_noise*randn(size(qp));
                        
%             const = 0.5 * norm(qpm, 'fro')^2;
            % find prediction error
            trueError_noise(ns_iter, idx, iter) = ...
                                evalSensorLocation(xsol_noise,Ad,eye(m),...
                                                xh, qpm, Qd, Rd);
            dmdError_noise(ns_iter, idx, iter) = ...
                                evalSensorLocation(xsolDMD_noise,Mdmd, Ldmd,...
                                            xh, qpm, Qdmd, Rd);
        end
    end
end
%%
trueColor = [0.8500 0.3250 0.0980];
dmdColor = [0.4940 0.1840 0.5560];

for ns_iter = 1:1%numel(nsList)
    figure('Position',[10 200 1000 800]);
    
    data1 = zeros(niter, numel(ntlong_list));
    for i = 1:numel(ntlong_list)
        data1(:,i) = trueError_noise(ns_iter,i,:);
    end
    
    data2 = zeros(niter, numel(ntlong_list));
    for i = 1:numel(ntlong_list)
        data2(:,i) = dmdError_noise(ns_iter,i,:);
    end

    b = boxplot(data1./data2, ...
            'Labels', {'0.005', '0.01', '0.02', '0.05'}, ...
            'BoxStyle', 'outline', ...
            'Colors', 'k', 'Symbol','w+');
    set(b,{'linew'},{2})
    
%     hold on
% 
%     b = boxplot(data2, ...
%             'Labels', {'100', '200', '300', '400'}, ...
%             'BoxStyle', 'outline', ...
%             'Colors', dmdColor);
% 	set(b,{'linew'},{2})    

    hold on

    jitter = 0.05;
    
    x1 = ones(length(data1),1);
    scatter(x1 + jitter*randn(size(data1(:,1))),data1(:,1)./data2(:,1), 20, ...
                'filled', 'k', 'o', 'MarkerEdgeColor', 'k');
    scatter(x1*2 + jitter*randn(size(data1(:,1))),data1(:,2)./data2(:,2), 20, 'filled', 'k', 'o', 'MarkerEdgeColor', 'k');
    scatter(x1*3 + jitter*randn(size(data1(:,1))),data1(:,3)./data2(:,3), 20, 'filled', 'k', 'o', 'MarkerEdgeColor', 'k');
    scatter(x1*4 + jitter*randn(size(data1(:,1))),data1(:,4)./data2(:,4), 20, 'filled', 'k', 'o', 'MarkerEdgeColor', 'k');

    hold on
    
    plot(0:5, ones(6,1), '-.k', 'linewidth', 1);
    
    hold off

    
%     box_vars = findall(gca,'Tag','Box');
%     leg = legend(box_vars([5,1]), {'full-order', 'DMD'});
%     set(leg,'Interpreter','latex','location','northwest');
    
    set(gca, 'FontSize', 30)
    %legend('true', 'dmd')
    xlabel('noise level','interpreter','latex');
    ylabel('relative estimation error','interpreter','latex');
    
end

h = gcf;
set(h,'PaperOrientation','landscape');                                                                                                     
set(h,'PaperUnits','normalized');                                                                                                          
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', 'figures/fig4_c.pdf')
% ylim([0.95, 1.1]);
% print(gcf, '-dpdf', 'figures/fig4_c_review.pdf')

%%
% noise level analysis
% horizon t = 400
nt = 401;
% prediction and reconstruction analysis
rng(0);
%nsList = [1, 2, 3, 4, 5];
nsList = [1];

%ntlong_list = [50, 100, 150, 200, 250, 300, 350, 400];
%ntlong_list = [100, 200, 300, 400];

noise_list = [0.005, 0.01, 0.02, 0.05];

niter = 20;

trueError_noise = zeros(numel(nsList), numel(noise_list), niter);
dmdError_noise = zeros(numel(nsList), numel(noise_list), niter);

for ns_iter = 1:numel(nsList)
    
    ns = nsList(ns_iter);
            
    xsol_noise = optimalSensors{ns_iter,1};
    xsolDMD_noise = optimalSensors{ns_iter,2};
    
    for idx = 1:numel(noise_list)

        curr_q_noise = noise_list(idx);
        curr_r_noise = noise_list(idx);
        Qd = curr_q_noise*(eye(m));

        % Q and R matrices
        Qdmd = Ldmd'*Qd*Ldmd;
        Rd = curr_r_noise*(eye(ns));        
        
        for iter = 1:niter

            disp([num2str(ns),'---',num2str(curr_q_noise),'---',num2str(iter)]);
            
            % simulate noisy data
            q0 = 1/5/(sqrt(2*pi))*exp(-(xh-8 + randn).^2/2/5^2);

            qp = zeros(size(q0,1),nt);
            qp(:,1) = q0;

            for i=2:nt
                    qp(:,i) = Ad*qp(:,i-1) + curr_q_noise*randn(size(m,1));
            end

            % measurement noise
            qpm = qp + curr_r_noise*randn(size(qp));
                        
%             const = 0.5 * norm(qpm, 'fro')^2;
            % find prediction error
            trueError_noise(ns_iter, idx, iter) = ...
                                evalSensorLocation(xsol_noise,Ad,eye(m),...
                                                xh, qpm, Qd, Rd);
            dmdError_noise(ns_iter, idx, iter) = ...
                                evalSensorLocation(xsolDMD_noise,Mdmd, Ldmd,...
                                            xh, qpm, Qdmd, Rd);
        end
    end
end
%%
trueColor = [0.8500 0.3250 0.0980];
dmdColor = [0.4940 0.1840 0.5560];

for ns_iter = 1:1%numel(nsList)
    figure('Position',[10 200 1000 800]);
    
    data1 = zeros(niter, numel(ntlong_list));
    for i = 1:numel(ntlong_list)
        data1(:,i) = trueError_noise(ns_iter,i,:);
    end
    
    data2 = zeros(niter, numel(ntlong_list));
    for i = 1:numel(ntlong_list)
        data2(:,i) = dmdError_noise(ns_iter,i,:);
    end

    b = boxplot(data1./data2, ...
            'Labels', {'0.005', '0.01', '0.02', '0.05'}, ...
            'BoxStyle', 'outline', ...
            'Colors', 'k', 'Symbol','w+');
    set(b,{'linew'},{2})
    
%     hold on
% 
%     b = boxplot(data2, ...
%             'Labels', {'100', '200', '300', '400'}, ...
%             'BoxStyle', 'outline', ...
%             'Colors', dmdColor);
% 	set(b,{'linew'},{2})

    hold on

    jitter = 0.05;
    
    x1 = ones(length(data1),1);
    scatter(x1 + jitter*randn(size(data1(:,1))),data1(:,1)./data2(:,1), 20, ...
                'filled', 'k', 'o', 'MarkerEdgeColor', 'k');
    scatter(x1*2 + jitter*randn(size(data1(:,1))),data1(:,2)./data2(:,2), 20, 'filled', 'k', 'o', 'MarkerEdgeColor', 'k');
    scatter(x1*3 + jitter*randn(size(data1(:,1))),data1(:,3)./data2(:,3), 20, 'filled', 'k', 'o', 'MarkerEdgeColor', 'k');
    scatter(x1*4 + jitter*randn(size(data1(:,1))),data1(:,4)./data2(:,4), 20, 'filled', 'k', 'o', 'MarkerEdgeColor', 'k');

    hold on
    
    plot(0:5, ones(6,1), '-.k', 'linewidth', 1);
    
    hold off

    
%     box_vars = findall(gca,'Tag','Box');
%     leg = legend(box_vars([5,1]), {'full-order', 'DMD'});
%     set(leg,'Interpreter','latex','location','northwest');
    
    set(gca, 'FontSize', 30)
    %legend('true', 'dmd')
    xlabel('noise level','interpreter','latex');
    ylabel('relative estimation error','interpreter','latex');
    
end

h = gcf;
set(h,'PaperOrientation','landscape');                                                                                                     
set(h,'PaperUnits','normalized');                                                                                                          
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', 'figures/fig4_d.pdf')

%%%%%%
%%
% noise level analysis
% horizon t = 200
nt = 201;
% prediction and reconstruction analysis
rng(0);
%nsList = [1, 2, 3, 4, 5];
nsList = [1];

%ntlong_list = [50, 100, 150, 200, 250, 300, 350, 400];
%ntlong_list = [100, 200, 300, 400];

noise_list = [0.005, 0.01, 0.02, 0.05];

niter = 20;

trueError_noise = zeros(numel(nsList), numel(noise_list), niter);
dmdError_noise = zeros(numel(nsList), numel(noise_list), niter);

for ns_iter = 1:numel(nsList)
    
    ns = nsList(ns_iter);
            
    xsol_noise = optimalSensors{ns_iter,1};
    xsolDMD_noise = optimalSensors{ns_iter,2};
    
    for idx = 1:numel(noise_list)

        curr_q_noise = noise_list(idx);
        curr_r_noise = noise_list(idx);
        Qd = 0.01*(eye(m));

        % Q and R matrices
        Qdmd = Ldmd'*Qd*Ldmd;
        Rd = curr_r_noise*(eye(ns));        
        
        for iter = 1:niter

            disp([num2str(ns),'---',num2str(curr_q_noise),'---',num2str(iter)]);
            
            % simulate noisy data
            q0 = 1/5/(sqrt(2*pi))*exp(-(xh-8 + randn).^2/2/5^2);

            qp = zeros(size(q0,1),nt);
            qp(:,1) = q0;

            for i=2:nt
                    qp(:,i) = Ad*qp(:,i-1) + curr_q_noise*randn(size(m,1));
            end

            % measurement noise
            qpm = qp + curr_r_noise*randn(size(qp));
            
%             const = 0.5 * norm(qpm, 'fro')^2;
            % find prediction error
            trueError_noise(ns_iter, idx, iter) = ...
                                evalSensorLocation(xsol_noise,Ad,eye(m),...
                                                xh, qpm, Qd, Rd);
            dmdError_noise(ns_iter, idx, iter) = ...
                                evalSensorLocation(xsolDMD_noise,Mdmd, Ldmd,...
                                            xh, qpm, Qdmd, Rd);
        end
    end
end
%%
trueColor = [0.8500 0.3250 0.0980];
dmdColor = [0.4940 0.1840 0.5560];

for ns_iter = 1:1%numel(nsList)
    figure('Position',[10 200 1000 800]);
    
    data1 = zeros(niter, numel(ntlong_list));
    for i = 1:numel(ntlong_list)
        data1(:,i) = trueError_noise(ns_iter,i,:);
    end
    
    data2 = zeros(niter, numel(ntlong_list));
    for i = 1:numel(ntlong_list)
        data2(:,i) = dmdError_noise(ns_iter,i,:);
    end

    b = boxplot(data1./data2, ...
            'Labels', {'0.005', '0.01', '0.02', '0.05'}, ...
            'BoxStyle', 'outline', ...
            'Colors', 'k', 'Symbol','w+');
    set(b,{'linew'},{2})
    
%     hold on
% 
%     b = boxplot(data2, ...
%             'Labels', {'100', '200', '300', '400'}, ...
%             'BoxStyle', 'outline', ...
%             'Colors', dmdColor);
% 	set(b,{'linew'},{2})    

    hold on

    jitter = 0.05;
    
    x1 = ones(length(data1),1);
    scatter(x1 + jitter*randn(size(data1(:,1))),data1(:,1)./data2(:,1), 20, ...
                'filled', 'k', 'o', 'MarkerEdgeColor', 'k');
    scatter(x1*2 + jitter*randn(size(data1(:,1))),data1(:,2)./data2(:,2), 20, 'filled', 'k', 'o', 'MarkerEdgeColor', 'k');
    scatter(x1*3 + jitter*randn(size(data1(:,1))),data1(:,3)./data2(:,3), 20, 'filled', 'k', 'o', 'MarkerEdgeColor', 'k');
    scatter(x1*4 + jitter*randn(size(data1(:,1))),data1(:,4)./data2(:,4), 20, 'filled', 'k', 'o', 'MarkerEdgeColor', 'k');

    hold on
    
    plot(0:5, ones(6,1), '-.k', 'linewidth', 1);

    hold off

    
%     box_vars = findall(gca,'Tag','Box');
%     leg = legend(box_vars([5,1]), {'full-order', 'DMD'});
%     set(leg,'Interpreter','latex','location','northwest');
    
    set(gca, 'FontSize', 30)
    %legend('true', 'dmd')
    xlabel('true noise level','interpreter','latex');
    ylabel('relative estimation error','interpreter','latex');
    
end

ylim([0.95, 1.1]);

h = gcf;
set(h,'PaperOrientation','landscape');                                                                                                     
set(h,'PaperUnits','normalized');                                                                                                          
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', 'figures/fig4_extra1.pdf')

%%
% noise level analysis
% horizon t = 400
nt = 401;
% prediction and reconstruction analysis
rng(0);
%nsList = [1, 2, 3, 4, 5];
nsList = [1];

%ntlong_list = [50, 100, 150, 200, 250, 300, 350, 400];
%ntlong_list = [100, 200, 300, 400];

noise_list = [0.005, 0.01, 0.02, 0.05];

niter = 20;

trueError_noise = zeros(numel(nsList), numel(noise_list), niter);
dmdError_noise = zeros(numel(nsList), numel(noise_list), niter);

for ns_iter = 1:numel(nsList)
    
    ns = nsList(ns_iter);
            
    xsol_noise = optimalSensors{ns_iter,1};
    xsolDMD_noise = optimalSensors{ns_iter,2};
    
    for idx = 1:numel(noise_list)

        curr_q_noise = noise_list(idx);
        curr_r_noise = noise_list(idx);
        Qd = 0.01*(eye(m));

        % Q and R matrices
        Qdmd = Ldmd'*Qd*Ldmd;
        Rd = curr_r_noise*(eye(ns));        
        
        for iter = 1:niter

            disp([num2str(ns),'---',num2str(curr_q_noise),'---',num2str(iter)]);
            
            % simulate noisy data
            q0 = 1/5/(sqrt(2*pi))*exp(-(xh-8 + randn).^2/2/5^2);

            qp = zeros(size(q0,1),nt);
            qp(:,1) = q0;

            for i=2:nt
                    qp(:,i) = Ad*qp(:,i-1) + curr_q_noise*randn(size(m,1));
            end

            % measurement noise
            qpm = qp + curr_r_noise*randn(size(qp));
                        
%             const = 0.5 * norm(qpm, 'fro')^2;
            % find prediction error
            trueError_noise(ns_iter, idx, iter) = ...
                                evalSensorLocation(xsol_noise,Ad,eye(m),...
                                                xh, qpm, Qd, Rd);
            dmdError_noise(ns_iter, idx, iter) = ...
                                evalSensorLocation(xsolDMD_noise,Mdmd, Ldmd,...
                                            xh, qpm, Qdmd, Rd);
        end
    end
end
%%
trueColor = [0.8500 0.3250 0.0980];
dmdColor = [0.4940 0.1840 0.5560];

for ns_iter = 1:1%numel(nsList)
    figure('Position',[10 200 1000 800]);
    
    data1 = zeros(niter, numel(ntlong_list));
    for i = 1:numel(ntlong_list)
        data1(:,i) = trueError_noise(ns_iter,i,:);
    end
    
    data2 = zeros(niter, numel(ntlong_list));
    for i = 1:numel(ntlong_list)
        data2(:,i) = dmdError_noise(ns_iter,i,:);
    end

    b = boxplot(data1./data2, ...
            'Labels', {'0.005', '0.01', '0.02', '0.05'}, ...
            'BoxStyle', 'outline', ...
            'Colors', 'k', 'Symbol','w+');
    set(b,{'linew'},{2})
    
%     hold on
% 
%     b = boxplot(data2, ...
%             'Labels', {'100', '200', '300', '400'}, ...
%             'BoxStyle', 'outline', ...
%             'Colors', dmdColor);
% 	set(b,{'linew'},{2})

    hold on

    jitter = 0.05;
    
    x1 = ones(length(data1),1);
    scatter(x1 + jitter*randn(size(data1(:,1))),data1(:,1)./data2(:,1), 20, ...
                'filled', 'k', 'o', 'MarkerEdgeColor', 'k');
    scatter(x1*2 + jitter*randn(size(data1(:,1))),data1(:,2)./data2(:,2), 20, 'filled', 'k', 'o', 'MarkerEdgeColor', 'k');
    scatter(x1*3 + jitter*randn(size(data1(:,1))),data1(:,3)./data2(:,3), 20, 'filled', 'k', 'o', 'MarkerEdgeColor', 'k');
    scatter(x1*4 + jitter*randn(size(data1(:,1))),data1(:,4)./data2(:,4), 20, 'filled', 'k', 'o', 'MarkerEdgeColor', 'k');

    hold on
    
    plot(0:5, ones(6,1), '-.k', 'linewidth', 1);
    
    hold off

    
%     box_vars = findall(gca,'Tag','Box');
%     leg = legend(box_vars([5,1]), {'full-order', 'DMD'});
%     set(leg,'Interpreter','latex','location','northwest');
    
    set(gca, 'FontSize', 30)
    %legend('true', 'dmd')
    xlabel('true noise level','interpreter','latex');
    ylabel('relative estimation error','interpreter','latex');
    
end

h = gcf;
set(h,'PaperOrientation','landscape');                                                                                                     
set(h,'PaperUnits','normalized');                                                                                                          
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', 'figures/fig4_extra2.pdf')


%%%%%%
%%
% rank analysis
% horizon t = 200
nt = 201;
% prediction and reconstruction analysis
rng(0);
%nsList = [1, 2, 3, 4, 5];
nsList = [1];

%ntlong_list = [50, 100, 150, 200, 250, 300, 350, 400];
%ntlong_list = [100, 200, 300, 400];

%noise_list = [0.005, 0.01, 0.02, 0.05];
noise = 0.01;
rank_list = [5, 10, 15, 20];

niter = 20;

trueError_rank = zeros(numel(nsList), numel(rank_list), niter);
dmdError_rank = zeros(numel(nsList), numel(rank_list), niter);

for ns_iter = 1:numel(nsList)
    
    ns = nsList(ns_iter);
            
    xsol_noise = optimalSensors{ns_iter,1};
    xsolDMD_noise = optimalSensors{ns_iter,2};

    for idx = 1:numel(rank_list)

        [Ldmd, Mdmd, ~] = get_low_rank(qm(:,2:end),qm(:,1:end-1),...
            rank_list(idx), 'dmd');
        
        curr_q_noise = noise;
        curr_r_noise = noise;
        Qd = curr_q_noise*(eye(m));

        % Q and R matrices
        Qdmd = Ldmd'*Qd*Ldmd;
        Rd = curr_r_noise*(eye(ns));        
        
        for iter = 1:niter

            disp([num2str(ns),'---',num2str(rank_list(idx)),'---',num2str(iter)]);
            
            % simulate noisy data
            q0 = 1/5/(sqrt(2*pi))*exp(-(xh-8 + randn).^2/2/5^2);

            qp = zeros(size(q0,1),nt);
            qp(:,1) = q0;

            for i=2:nt
                    qp(:,i) = Ad*qp(:,i-1) + curr_q_noise*randn(size(m,1));
            end

            % measurement noise
            qpm = qp + curr_r_noise*randn(size(qp));
                        
%             const = 0.5 * norm(qpm, 'fro')^2;
            % find prediction error
            trueError_rank(ns_iter, idx, iter) = ...
                                evalSensorLocation(xsol_noise,Ad,eye(m),...
                                                xh, qpm, Qd, Rd);
            dmdError_rank(ns_iter, idx, iter) = ...
                                evalSensorLocation(xsolDMD_noise,Mdmd, Ldmd,...
                                            xh, qpm, Qdmd, Rd);
        end
    end
end
%%
trueColor = [0.8500 0.3250 0.0980];
dmdColor = [0.4940 0.1840 0.5560];

for ns_iter = 1:1%numel(nsList)
    figure('Position',[10 200 1000 800]);
    
    data1 = zeros(niter, numel(rank_list));
    for i = 1:numel(rank_list)
        data1(:,i) = trueError_rank(ns_iter,i,:);
    end
    
    data2 = zeros(niter, numel(rank_list));
    for i = 1:numel(rank_list)
        data2(:,i) = dmdError_rank(ns_iter,i,:);
    end

    b = boxplot(data1./data2, ...
            'Labels', {'5', '10', '15', '20'}, ...
            'BoxStyle', 'outline', ...
            'Colors', 'k', 'Symbol','w+');
    set(b,{'linew'},{2})
    
%     hold on
% 
%     b = boxplot(data2, ...
%             'Labels', {'100', '200', '300', '400'}, ...
%             'BoxStyle', 'outline', ...
%             'Colors', dmdColor);
% 	set(b,{'linew'},{2})    

    hold on

    jitter = 0.05;
    
    x1 = ones(length(data1),1);
    scatter(x1 + jitter*randn(size(data1(:,1))),data1(:,1)./data2(:,1), 20, ...
                'filled', 'k', 'o', 'MarkerEdgeColor', 'k');
    scatter(x1*2 + jitter*randn(size(data1(:,1))),data1(:,2)./data2(:,2), 20, 'filled', 'k', 'o', 'MarkerEdgeColor', 'k');
    scatter(x1*3 + jitter*randn(size(data1(:,1))),data1(:,3)./data2(:,3), 20, 'filled', 'k', 'o', 'MarkerEdgeColor', 'k');
    scatter(x1*4 + jitter*randn(size(data1(:,1))),data1(:,4)./data2(:,4), 20, 'filled', 'k', 'o', 'MarkerEdgeColor', 'k');
    %scatter(x1*4 + jitter*randn(size(data1(:,1))),ones(size(data1(:,1))), 20, 'filled', 'k', 'o', 'MarkerEdgeColor', 'k');

    hold on
    
    plot(0:5, ones(6,1), '-.k', 'linewidth', 1);
    
    hold off

    
%     box_vars = findall(gca,'Tag','Box');
%     leg = legend(box_vars([5,1]), {'full-order', 'DMD'});
%     set(leg,'Interpreter','latex','location','northwest');
    
    set(gca, 'FontSize', 30);
    %legend('true', 'dmd')
    xlabel('rank','interpreter','latex');
    ylabel('relative estimation error','interpreter','latex');
    
end

h = gcf;
set(h,'PaperOrientation','landscape');                                                                                                     
set(h,'PaperUnits','normalized');                                                                                                          
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', 'figures/fig4_e.pdf')

%%
% noise level analysis
% horizon t = 400
nt = 401;
% prediction and reconstruction analysis
rng(0);
%nsList = [1, 2, 3, 4, 5];
nsList = [1];

%ntlong_list = [50, 100, 150, 200, 250, 300, 350, 400];
%ntlong_list = [100, 200, 300, 400];

%noise_list = [0.005, 0.01, 0.02, 0.05];
noise = 0.01;
rank_list = [5, 10, 15, 20];

niter = 20;

trueError_rank = zeros(numel(nsList), numel(rank_list), niter);
dmdError_rank = zeros(numel(nsList), numel(rank_list), niter);

for ns_iter = 1:numel(nsList)
    
    ns = nsList(ns_iter);
            
    xsol_noise = optimalSensors{ns_iter,1};
    xsolDMD_noise = optimalSensors{ns_iter,2};

    for idx = 1:numel(rank_list)

        [Ldmd, Mdmd, ~] = get_low_rank(qm(:,2:end),qm(:,1:end-1),...
            rank_list(idx), 'dmd');
        
        curr_q_noise = noise;
        curr_r_noise = noise;
        Qd = curr_q_noise*(eye(m));

        % Q and R matrices
        Qdmd = Ldmd'*Qd*Ldmd;
        Rd = curr_r_noise*(eye(ns));        
        
        for iter = 1:niter

            disp([num2str(ns),'---',num2str(rank_list(idx)),'---',num2str(iter)]);
            
            % simulate noisy data
            q0 = 1/5/(sqrt(2*pi))*exp(-(xh-8 + randn).^2/2/5^2);

            qp = zeros(size(q0,1),nt);
            qp(:,1) = q0;

            for i=2:nt
                    qp(:,i) = Ad*qp(:,i-1) + curr_q_noise*randn(size(m,1));
            end

            % measurement noise
            qpm = qp + curr_r_noise*randn(size(qp));
                        
%             const = 0.5 * norm(qpm, 'fro')^2;
            % find prediction error
            trueError_rank(ns_iter, idx, iter) = ...
                                evalSensorLocation(xsol_noise,Ad,eye(m),...
                                                xh, qpm, Qd, Rd);
            dmdError_rank(ns_iter, idx, iter) = ...
                                evalSensorLocation(xsolDMD_noise,Mdmd, Ldmd,...
                                            xh, qpm, Qdmd, Rd);
        end
    end
end
%%
trueColor = [0.8500 0.3250 0.0980];
dmdColor = [0.4940 0.1840 0.5560];

for ns_iter = 1:1%numel(nsList)
    figure('Position',[10 200 1000 800]);
    
    data1 = zeros(niter, numel(rank_list));
    for i = 1:numel(rank_list)
        data1(:,i) = trueError_rank(ns_iter,i,:);
    end
    
    data2 = zeros(niter, numel(rank_list));
    for i = 1:numel(rank_list)
        data2(:,i) = dmdError_rank(ns_iter,i,:);
    end

    b = boxplot(data1./data2, ...
            'Labels', {'5', '10', '15', '20'}, ...
            'BoxStyle', 'outline', ...
            'Colors', 'k', 'Symbol','w+');
    set(b,{'linew'},{2})
    
%     hold on
% 
%     b = boxplot(data2, ...
%             'Labels', {'100', '200', '300', '400'}, ...
%             'BoxStyle', 'outline', ...
%             'Colors', dmdColor);
% 	set(b,{'linew'},{2})    

    hold on

    jitter = 0.05;
    
    x1 = ones(length(data1),1);
    scatter(x1 + jitter*randn(size(data1(:,1))),data1(:,1)./data2(:,1), 20, ...
                'filled', 'k', 'o', 'MarkerEdgeColor', 'k');
    scatter(x1*2 + jitter*randn(size(data1(:,1))),data1(:,2)./data2(:,2), 20, 'filled', 'k', 'o', 'MarkerEdgeColor', 'k');
    scatter(x1*3 + jitter*randn(size(data1(:,1))),data1(:,3)./data2(:,3), 20, 'filled', 'k', 'o', 'MarkerEdgeColor', 'k');
    scatter(x1*4 + jitter*randn(size(data1(:,1))),data1(:,4)./data2(:,4), 20, 'filled', 'k', 'o', 'MarkerEdgeColor', 'k');
    %scatter(x1*4 + jitter*randn(size(data1(:,1))),ones(size(data1(:,1))), 20, 'filled', 'k', 'o', 'MarkerEdgeColor', 'k');

    hold on
    
    plot(0:5, ones(6,1), '-.k', 'linewidth', 1);
    
    hold off
    
%     box_vars = findall(gca,'Tag','Box');
%     leg = legend(box_vars([5,1]), {'full-order', 'DMD'});
%     set(leg,'Interpreter','latex','location','northwest');
    
    set(gca, 'FontSize', 30);
    %legend('true', 'dmd')
    xlabel('rank','interpreter','latex');
    ylabel('relative estimation error','interpreter','latex');
    
end

h = gcf;
set(h,'PaperOrientation','landscape');                                                                                                     
set(h,'PaperUnits','normalized');                                                                                                          
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', 'figures/fig4_f.pdf')