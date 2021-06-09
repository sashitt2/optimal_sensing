function [C, dCdxs] = getCfromSensorloc(xs, xh, sigma)

%     C = exp(-(xh-xs(i)).^2/2/sigma^2)'*M;
%     dCdxs = C.*(xh - xs(i))'/sigma^2;

if numel(sigma) == 1
    sigma = sigma * ones(size(xs));
end

M = diag(([diff(xh);0]+[0;diff(xh)])/2);

C = zeros(numel(xs),numel(xh));
dCdxs = zeros(numel(xs),numel(xh));

for i = 1:size(xs,1)
    C(i,:) = exp(-(xh-xs(i)).^2/2/sigma(i)^2)'*M;
    dCdxs(i,:) = C(i,:).*(xh - xs(i))'/sigma(i)^2;
end

end