function error = evalSSPORsensor(r,ns,q)

[psi, ~, ~] = svd(q, 'econ');

[~,~,pivot] = qr(psi(:,1:r)','vector');
pivot = pivot(1:ns);

qqr = psi(:,1:r)*(psi(pivot,1:r)\q(pivot,:));

error = 0.5 * norm(qqr - q, 'fro')^2;

end