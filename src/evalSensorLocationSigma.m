function objValue = evalSensorLocationSigma(xs, sigma,A0,L0,xh,q,Q,R)

% data series
q0 = q(:,1);
q = q(:,2:end);

% observation matrix
%C = exp(-(xh-xs0).^2/2/sigma^2)'*M;
C = getCfromSensorlocwidth(xs, xh, sigma);

objValue = getObjInfinite(A0, q0, q, Q, R, C, L0);

end