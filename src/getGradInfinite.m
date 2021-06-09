function [gradX, gradA, gradC, gradL] = getGradInfinite(A,x0,x,Q,R,C,L)


r = size(A,1);
n = size(x,2);

CL = C*L;

%% forward Kalman filter

x_sub0 = L'*x0;
x_sub = zeros(r,1);
x_est = zeros(r,1);

assert(size(Q,1) == size(Q,2), 'Q matrix is not square');

if size(Q,1) == 1
    LtQL = Q*eye(r);
elseif size(Q,1) == r
    LtQL = Q;
    Q = Q(1,1);
else
    disp(['size of Q = ', num2str(size(Q,1)), ' not consistent']);
    return
end

% offline computation

Sigma = dare(A',CL',LtQL,R);

S = CL*Sigma*CL' + R;
invS = inv(S);

K = Sigma*CL'/S;

% online computation
for i=0:1:n-1
    
    if(i==0)
        x_est(:,i+1) = A*x_sub0;
    else
        x_est(:,i+1) = A*x_sub(:,i);
    end
    
    x_sub(:,i+1) = x_est(:,i+1) + K*(C*(x(:,i+1) - L*x_est(:,i+1)));

end

%r_est = (x - x_sub);
r_est = (L'*x - x_sub);

%% adjoint kalman filter

eta_est = zeros(size(x_est));
eta_sub = zeros(size(x_sub));
H = zeros(size(K));

%% offline computation
eta_sub(:,n) = -r_est(:,n);

for i=n:-1:1

    eta_est(:,i) = eta_sub(:,i) - CL'*K'*eta_sub(:,i);

    if i > 1
        eta_sub(:,i-1) = A'*eta_est(:,i) - r_est(:,i-1);
    end
    
    %H = H + eta_sub(:,i)*(x(:,i) - x_est(:,i))'*CL';
    H = H + eta_sub(:,i)*(x(:,i) - L*x_est(:,i))'*C';
end

eta_sub0 = A'*eta_est(:,1);

%% online computation
SitC = S'\CL;

% tic
% rhs = (CL' * SitC * Sigma' - eye(r)) * H * SitC;
% 
% lhsMat = (-kron(eye(r),eye(r)) + kron(transp(A),A') ...
%           -kron(transp(A*Sigma'*CL'*SitC),A') ...
%           -kron(transp(A),CL' * SitC * Sigma' * A') ...
%           +kron(transp(A*Sigma'*CL'*SitC),CL' * SitC * Sigma' * A'));
% 
% LambdaVec = lhsMat\rhs(:);
% 
% Lambda = reshape(LambdaVec, size(Sigma));

Z = (eye(r) - CL'*K');
Lambda = dlyap(Z*A', Z*H*SitC);

%W = SitC * Sigma'* (A' * Lambda * A * Sigma' * CL' - H) *invS';
W = K' * (A' * Lambda * A * Sigma' * CL' - H) * invS';

%% gradient computation

gradA = eta_est(:,1)*x_sub0';

for i=n-1:-1:1
    gradA = gradA + eta_est(:,i+1)*x_sub(:,i)';
end

gradA = gradA + Lambda*A*Sigma' * ( eye(r) - CL'*SitC*Sigma' ) ...
              + Lambda'*A*Sigma * ( eye(r) - SitC'*CL*Sigma  );

gradL = 0;
for i=n:-1:1
    %gradL = gradL + C'*K'*eta_sub(:,i)*(x(:,i) - x_est(:,i))';
    gradL = gradL - C'*K'*eta_sub(:,i)*x_est(:,i)';
end

gradL = gradL + C'/S*H'*Sigma + C'*W*CL*Sigma' + C'*W'*CL*Sigma ...
        + Q*L*Lambda' + Q'*L*Lambda ...
        -(  C'*SitC*Sigma*A'*Lambda'*A*Sigma ...
          + C'*SitC*Sigma'*A'*Lambda*A*Sigma');

gradL = gradL + x0*eta_sub0';      

for i = 1:n
    gradL = gradL - x(:,i)*x_sub(:,i)';
end

for i = 1:n
    gradL = gradL + (L*x_sub(:,i))*x_sub(:,i)';
end


gradC = 0;
for i=n:-1:1
    %gradC = gradC + K'*eta_sub(:,i)*(x(:,i) - x_est(:,i))'*L';
    gradC = gradC + K'*eta_sub(:,i)*(x(:,i) - L * x_est(:,i))';
end

gradC = gradC + S\H'*Sigma*L' + W*CL*Sigma'*L' + W'*CL*Sigma*L' ...
        - ( SitC*Sigma*A'*Lambda'*A*Sigma*L' + ...
            SitC*Sigma'*A'*Lambda*A*Sigma'*L' );


gradX = L*eta_sub0;

end
