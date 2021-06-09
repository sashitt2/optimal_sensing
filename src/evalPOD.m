function error = evalPOD(r, q)

[psi, ~, ~] = svd(q, 'econ');

qpod = psi(:,1:r)*(psi(:,1:r)'*q);

error = 0.5 * norm(qpod - q, 'fro')^2;

end