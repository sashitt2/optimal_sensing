function error = evalSSPORsensorloc(r, xs, xh, q)

[psi, ~, ~] = svd(q, 'econ');

% [~,~,pivot] = qr(psi(:,1:r)','vector');
% pivot = pivot(1:ns);

pivot = zeros(size(xs));
for i = 1:size(xs,1)
    [~, minloc] = min(abs(xh - xs(i)));
    pivot(i) = minloc;
end

qqr = psi(:,1:r)*(psi(pivot,1:r)\q(pivot,:));

error = 0.5 * norm(qqr - q, 'fro')^2;

end