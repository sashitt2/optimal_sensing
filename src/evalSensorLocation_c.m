function objValue = evalSensorLocation_c(xs,A0,L0,B,u,xh,q,Q,R,sigma)

% mass matrix and sigma (sensor matrix)
if nargin < 10
    sigma = 2;
end

% data series
q0 = q(:,1);
q = q(:,2:end);

% observation matrix
%C = exp(-(xh-xs0).^2/2/sigma^2)'*M;
C = getCfromSensorloc(xs, xh, sigma);

objValue = getObjInfinite_c(A0, q0, q, Q, R, C, L0, B, u);

end