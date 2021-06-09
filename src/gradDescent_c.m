function [xs, A, L, objList] = gradDescent_c(xs0,A0,L0,B,u,xh,q,Q,R, ...
                                update,maxIter,step_init,sigma)

if nargin < 13
    sigma = 2;
end

if nargin < 12
    step_init = 1;
end

if nargin < 11
    %maxIter = 1000;
    maxIter = 25;
end

if nargin < 10
    update = [1, 1, 1];
end

if update(1)
    disp('updating xs');
end
if update(2)
    disp('updating A');
end
if update(3)
    disp('updating L');
end

% data series
q0 = q(:,1);
q = q(:,2:end);

% initialize initial guess
A = A0;
L = L0;
xs = xs0;

%C = exp(-(xh-xs).^2/2/sigma^2)'*M;
[C, ~] = getCfromSensorloc(xs, xh, sigma);

% initialize previous iteration
Alast = A0;
Llast = L0;
xslast = xs0;
Clast = C;

valueOld = getObjInfinite_c(A, q0, q, Q, R, C, L, B, u);
disp(['objective value of initial guess is ', num2str(valueOld)]);
%%

%print the reporting header
printInfo();

%get an initial step-size
%maxIter = 50;
%step_init = 1;
step = step_init;
% have to implement armijo rule here

%objList = zeros(maxIter,1);
objList = [];

%intialize convergence test variables
iterCount = 0;
relError  = 1;

tol = 10^(-4);
%%
while( abs(relError) > tol && iterCount < maxIter)

    %[~, dgdA, dgdC, dgdL] = getGradInfinite(A, q0, q, Q, R, C, L);
    [~, dgdA, dgdC, dgdL] = getGradInfinite_c(A, q0, q, Q, R, C, L, B, u);
    
    % update xs and C
    if update(1)

        %dCdxs = C.*(xh - xs)'/sigma^2;
        [~, dCdxs] = getCfromSensorloc(xs, xh, sigma);
                
        %dgdxs = sum(conj(dCdxs).*dgdC, 2);
        dgdxs = sum(conj(dCdxs).*dgdC, 2);

        xs = xs - step*real(dgdxs);
        %C = exp(-(xh-xs).^2/2/sigma^2)'*M;
        [C, ~] = getCfromSensorloc(xs, xh, sigma);
    end
    
    % update A
    if update(2)
        A = A - step*dgdA;
    end
    
    % update L
    if update(3)
        Ltemp = L - step*(dgdL - L*dgdL'*L);
        
        % retraction for L (Steifel manifold)
        [Qtemp,Rtemp] = qr(Ltemp,0);
        L = Qtemp * diag(sign(sign(diag(Rtemp))+.5));
    end
            
    valueNew = getObjInfinite_c(A, q0, q, Q, R, C, L, B, u);    
    
    %calculate relative improvement
    relError = (valueOld-valueNew)/abs(valueOld);
    
    %print reporting information
    if update(1)
        printInfo(iterCount,valueNew,norm(dgdxs),relError,step);
    else
        printInfo(iterCount,valueNew,norm(dgdA),relError,step);
    end
    
    if(relError <= 0 || step == 0 || isnan(valueNew))
        %not improving.  Force stop.
        relError = 1;
        
        % go back to previous iteration
        A = Alast; 
        L = Llast;
        
        xs = xslast;
        C = Clast;
        
        % decrease step size by factor of 2 (Armijo)
        step = step/2;
    else
        % move forward
        Alast = A;
        Llast = L;
        
        xslast = xs;
        Clast = C;
        
        % re-initialize the step size for next iteration
        step = step_init;
        
        % record the value
        %objList(iterCount + 1) = valueNew;
        objList = [objList; valueNew];
        
        %update counter
        iterCount = iterCount + 1;
        valueOld = valueNew; %for next pass
    end

end

end