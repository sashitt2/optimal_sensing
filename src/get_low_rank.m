function [L,M,R] = get_low_rank(A,B,r,method)
% X is the low rank approximation
% A,B are the data matrices
% A = X B
% r is the rank of the approximation
% method says which method to use

switch method
    
    case 'svd'
        % here A and B are both the true dynamical matrix
        disp(['performing svd']);
        [U,S,V] = svd(A,'econ');
        %X = U(:,1:r)*S(1:r,1:r)*V(:,1:r)';
        
    case 'dmd'
        % here A and B are the data matrices
        disp(['performing dmd']);
        [U,S,V] = svd(B,'econ');
        Ut = U(:,1:r);
        Vt = V(:,1:r);
        St = S(1:r,1:r);
        M = (Ut'*A)*Vt*inv(St);
        %[L,M] = omd(A,B,r,Ut,'dmd');
        %X = L*M*L';
        L = Ut;
        R = L;
        
    case 'tdmd'
        % here A and B are the data matrices
        disp(['performing tdmd']);
        Z = [B; A];
        [~, ~,Vf] = svd(Z, 'econ');
        A_bar = (A * Vf) * Vf';
        B_bar = (B * Vf) * Vf';
        norm(A - A_bar, 'fro')
        norm(B - B_bar, 'fro')
        [U_bar,S_bar,V_bar] = svd(B_bar,'econ');
        Ut = U_bar(:,1:r);
        Vt = V_bar(:,1:r);
        St = S_bar(1:r,1:r);
        M = (Ut'*A_bar)*Vt*inv(St);
        %[L,M] = omd(A,B,r,Ut,'dmd');
        %X = L*M*L';
        L = Ut;
        R = L;        
        
    case 'omd'
        % here A and B are the data matrices
        disp(['performing omd']);
        [U,S,V] = svd(B,'econ');
        Ut = U(:,1:r);
        %[L,M] = omd(A,B,r,Ut,'alternating');
        %[L,M] = omd(A,B,r,L,'conjgrad');
        [L,M,L] = lmlOpt_trustregion(A,B,r,Ut);
        %X = L*M*L';
        R = L;
        
    case 'ldr'
        % here A and B are the data matrices
        disp(['performing ldr']);
        [L,D,R] = ldrOpt_subProject(A,B,r);
        %X = L*D*R';
        M = D;
        
    case 'hyb'
        % here A and B are the data matrices
        disp(['performing hyb']);
        [L,D,R] = ldrOpt_hybrid(A,B,r);
        %X = L*D*R';
        M = D;        
        
    case 'olr'
        % here A and B are the data matrices
        disp(['performing olr']);
        [L,D,R] = ldrOpt_subProject(A,B,r);
        [L,D,R] = ldrOpt_gradient(A,B,L,R);
        %X = L*D*R';
        M = D;
        
    case 'conj'
        % here A and B are the data matrices
        disp(['performing ldr with conj']);
        [L,D,R] = ldrOpt_subProject(A,B,r);
        x0.L = L;
        x0.R = R;
        %[L,D,R] = ldrOpt_conjgrad(A,B,r,x0);
        [L,D,R] = ldrOpt_trustregion(A,B,r,x0);
        %X = L*D*R';
        M = D;

end

end