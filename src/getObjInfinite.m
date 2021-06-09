function J = getObjInfinite(A,x0,x,Q,R,C,L)

r = size(A,1);
n = size(x,2);

CL = C*L;

%% Kalman filter

x_sub0 = L'*x0;
x_sub = zeros(r,n);
x_est = zeros(r,n);

%%
% offline computation

try
    P= dare(A',CL',Q,R);
catch
    warning('Observer not feasible');
    J = realmax;
    return
end

%disp(['error is ', num2str(norm(P-A*P*A'-Q+A*P*C'*inv(C*P*C'+R)*C*P*A'))])

K = P*CL'/(CL*P*CL' + R);

% Sigma = dare(A',C',Q,R);
% 
% S = C*Sigma*C' + R;
%
% K = Sigma*C'/S;
%%
% online computation
for i=0:1:n-1

    if(i==0)
        x_est(:,i+1) = A*x_sub0;
    else
        x_est(:,i+1) = A*x_sub(:,i);
    end
    
    %x_sub(:,i+1) = x_est(:,i+1) + K*CL*(x(:,i+1) - x_est(:,i+1));
    x_sub(:,i+1) = x_est(:,i+1) + K*(C*(x(:,i+1) - L*x_est(:,i+1)));
end

%r_est = (x - x_sub);
r_est = (x - L*x_sub);

%%
% compute objective function

J = 0;

for i=1:n
     J = J + 0.5*(r_est(:,i)'*r_est(:,i));
%      J = J + 0.5*(x(:,i)'*x(:,i) + x_sub(:,i)'*x_sub(:,i) ...
%                     - 2 * x_sub(:,i)' * (L' * x(:,i)));
end

end
