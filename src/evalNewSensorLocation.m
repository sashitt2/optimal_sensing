function objValue = evalNewSensorLocation(xs,A0,L0,xh,q,Q,R)

% data series
q0 = q(:,1);
q = q(:,2:end);

% observation matrix
%C = exp(-(xh-xs0).^2/2/sigma^2)'*M;
%C = getCfromSensorloc(xs, xh, sigma);

C = zeros(size(xs,1), size(q,1));
for i = 1:size(xs,1)
    [~, minloc] = min(abs(xh - xs(i)));
    C(i, minloc) = 1;
end

objValue = getObjInfinite(A0, q0, q, Q, R, C, L0);

end