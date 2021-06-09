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
podColor = [0.12, 0.467, 0.706];
ssporBasicColor = [0.84, 0.153, 0.157];
ssporColor = [0.17254901960784313, 0.6274509803921569, 0.1725490196078431];
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
%%
% data has been generated. Now we look at the optimal sensor location for
% varying number of sensor locations

nsList = [1, 2, 3, 4, 5];
sigmaList = [1, 0.5, 0.25];
stepSizeList = [1, 0.8, 0.6];
numIterList = [50, 30, 10];
rng(1);

optimalSensors = cell(numel(nsList), 3);
timeTaken = zeros(numel(nsList), 3);
reconError = cell(numel(nsList),3);

Qd = q_noise*(eye(m));

r = 10;
[Ldmd, Mdmd, ~] = get_low_rank(qm(:,2:end),qm(:,1:end-1),r,'dmd');

update = [1, 0 , 0];

for iter = 1:numel(nsList)

    ns = nsList(iter);
    
    xs0 = iter * randn(ns, 1);
    
    Qdmd = Ldmd'*Qd*Ldmd;
    Rd = r_noise*(eye(ns));

    disp(['DMD dynamics optimal sensor location -> ', num2str(iter)]);
    tic
    curr_sigma = sigmaList(1);
    while true
        [xsolDMD, ~, ~, ~] = gradDescent(xs0, Mdmd, Ldmd, ...
                               xh, qm, Qdmd, Rd, update, 50, ...
                                        curr_sigma, curr_sigma);
        if curr_sigma <= 0.25
            break
        else
            xs0 = xsolDMD;        
            curr_sigma = 0.9 * curr_sigma;
        end
    end
    timeDMD = toc;
    
    disp(['SSPOR optimal sensor location -> ', num2str(iter)]);    
    tic
    xsolSSPOR = getSSPOR(xh, r, ns, qm);
    timeSSPOR = toc;
    
    disp(['SSPOR+DMD optimal sensor location -> ', num2str(iter)]);    
    tic
    xs0 = getSSPOR(xh, r, ns, qm);
    [xsolBOTH, ~, ~, ~] = gradDescent(xs0, Mdmd, Ldmd, ...
                            xh, qm, Qdmd, Rd, update, numIterList(end), ...
                            stepSizeList(end-1), sigmaList(end));    
    timeBOTH = toc;
    
    timeTaken(iter,1) = timeDMD;
    timeTaken(iter,2) = timeSSPOR;
    timeTaken(iter,3) = timeBOTH;
    
    optimalSensors{iter,1} = xsolDMD;
    optimalSensors{iter,2} = xsolSSPOR;
    optimalSensors{iter,3} = xsolBOTH;
    
    reconError{iter,1} = evalNewSensorLocation(xsolDMD, Mdmd, Ldmd,...
                                        xh, qm, Qdmd, Rd);
    reconError{iter,2} = evalNewSensorLocation(xsolSSPOR, Mdmd, Ldmd,...
                                        xh, qm, Qdmd, Rd);                                    
    reconError{iter,3} = evalNewSensorLocation(xsolBOTH, Mdmd, Ldmd,...
                                        xh, qm, Qdmd, Rd);                                    

end
%%
sspor_error = zeros(numel(nsList),1);
pod_error = zeros(numel(nsList),1);
for i = 1: numel(nsList)
    sspor_error(i) = evalSSPORsensor(r, nsList(i), qm);
    pod_error(i) = evalPOD(nsList(i), qm);
end
%%
% random sensor placement results
nIter = 20;
randomSensorError = zeros(nIter, numel(nsList));

rng(5);
for i = 1:numel(nsList)
    for iter = 1:nIter
        % randomly sample the amplification region
        xs = -8.6 + 8.6*(rand(i,1));
        
        Rd = r_noise*(eye(i));
        
        randomSensorError(iter,i) = evalNewSensorLocation(xs, Mdmd, Ldmd,...
                                        xh, qm, Qdmd, Rd);        

        
        %randomSensorError(iter, i) = evalSSPORsensorloc(r, xs, xh, qm);

    end
    disp(['finished ',num2str(i)]);
end

%% PLOTS
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
                        'MarkerEdgeColor', dmdColor);
    end
    
    dmdSensors = optimalSensors{i,2};
    for j = 1:numel(trueSensors)
        plot(dmdSensors(j), yvalue, 'x', 'MarkerSize', 20, 'LineWidth', 3, ...
                        'MarkerEdgeColor', ssporBasicColor);
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
%leg = legend(line_vars([36,32,29]), {'KF-DMD', 'SSPOR', 'Hybrid'});
leg = legend(line_vars([31,30]), {'KF-DMD', 'SSPOR'});
set(leg,'Interpreter','latex','location','northeast');
box on

h = gcf;
set(h,'PaperOrientation','landscape');                                                                                                     
set(h,'PaperUnits','normalized');                                                                                                          
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', 'figures/fig2_a.pdf')
%%
nsList = [1, 2, 3, 4, 5];
figure('Position',[10 200 1000 800]);

semilogy(nsList, timeTaken(:,1), '-ok', 'MarkerSize', 15, 'LineWidth', 3, ...
                'Color', dmdColor)
hold on
semilogy(nsList, timeTaken(:,2) ,'-xr', 'MarkerSize', 15, 'LineWidth', 3, ...
                'Color', ssporBasicColor)
set(gca, 'FontSize', 30)
leg = legend('KF-DMD', 'SSPOR');
set(leg,'Interpreter','latex','location','northeast');
xlabel('number of sensors','interpreter','latex');
ylabel('time taken (s)','interpreter','latex');

h = gcf;
set(h,'PaperOrientation','landscape');                                                                                                     
set(h,'PaperUnits','normalized');                                                                                                          
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', 'figures/fig2_c.pdf')

%%

noise_list = [0.005, 0.01, 0.02, 0.05, 0.1];

ssporError_noise = zeros(numel(noise_list),1);
dmdError_noise = zeros(numel(noise_list),1);

ns = 1;
for idx = 1:numel(noise_list)

    rng(0);
    q = zeros(size(q0,1),nt);
    q(:,1) = q0;

    q_noise = noise_list(idx);
    for i=2:nt
            q(:,i) = Ad*q(:,i-1) + q_noise*randn(size(m,1));
    end

    % measurement noise
    qm = q + r_noise*randn(size(q));

    [Ldmd, Mdmd, ~] = get_low_rank(qm(:,2:end),qm(:,1:end-1),r,'dmd');

    xs0 = iter * randn(1, 1);
    
    Qdmd = Ldmd'*Qd*Ldmd;
    Rd = r_noise*(eye(1));
    
    curr_sigma = sigmaList(1);
    while true
        [xsolDMD, ~, ~, ~] = gradDescent(xs0, Mdmd, Ldmd, ...
                               xh, qm, Qdmd, Rd, update, 50, ...
                                        curr_sigma, curr_sigma);
        if curr_sigma <= 0.25
            break
        else
            xs0 = xsolDMD;        
            curr_sigma = 0.9 * curr_sigma;
        end
    end

    xsolSSPOR = getSSPOR(xh, r, 1, qm);
    
    const = 0.5 * norm(qm, 'fro')^2;
    
    dmdError_noise(idx,1) = evalNewSensorLocation(xsolDMD, Mdmd, Ldmd,...
                                        xh, qm, Qdmd, Rd) / const;
    ssporError_noise(idx,1) = evalSSPORsensor(r, 1, qm) / const;    
end
%%

figure('Position',[10 200 1000 800]);

semilogy(noise_list, dmdError_noise, '-o', 'MarkerSize', 15, 'LineWidth', 3, ...
                'Color', dmdColor);
hold on
semilogy(noise_list, ssporError_noise, '-x', 'MarkerSize', 15, 'LineWidth', 3, ...
                'Color', ssporBasicColor);

set(gca, 'FontSize', 30)
xlabel('system noise level','interpreter','latex');
ylabel('normalized estimation error','interpreter','latex');

leg = legend('KF-DMD', 'SSPOR');
set(leg,'Interpreter','latex','location','northeast');

h = gcf;
set(h,'PaperOrientation','landscape');                                                                                                     
set(h,'PaperUnits','normalized');                                                                                                          
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', 'figures/fig2_d.pdf')

%%
figure('Position',[10 200 1000 800]);

reconErrorMat = cell2mat(reconError);

const = 0.5 * norm(qm, 'fro')^2;

plot(nsList, reconErrorMat(:,1) / const, '-o', 'MarkerSize', 15, 'LineWidth', 3, ...
                'Color', dmdColor)
hold on

plot(nsList, sspor_error / const, '-x', 'MarkerSize', 15, 'LineWidth', 3, ...
                'Color', ssporBasicColor)
hold on
plot(nsList, pod_error / const ,'-^', 'MarkerSize', 15, 'LineWidth', 3, ...
                'Color', podColor)            
            
set(gca, 'FontSize', 30)

hold on

data = randomSensorError / const;

line_vars = findall(gca, 'Type', 'Line');
leg = legend(line_vars([1,2,3]), {'POD', 'SSPOR', 'KF-DMD'});
set(leg,'Interpreter','latex','location','northeast');

xlabel('number of sensors','interpreter','latex');
ylabel('normalized estimation error','interpreter','latex');

h = gcf;
set(h,'PaperOrientation','landscape');                                                                                                     
set(h,'PaperUnits','normalized');                                                                                                          
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', 'figures/fig2_b.pdf')