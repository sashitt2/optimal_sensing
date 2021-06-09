function xsensor = getSSPOR(xh, r, p, q)

[psi, ~, ~] = svd(q, 'econ');

[~,~,pivot] = qr(psi(:,1:r)','vector');
pivot = pivot(1:p);

xsensor = xh(pivot);

end