function [ Q, zEnd ] = pr_GaussNewton( pr )
% Phase retrieval using the Gauss-Newton algorithm

    z = pr.z0;
    J = zeros(size(pr.A)) ;
    
    Q(1) = abs(z.'*conj(pr.trg.x))^2/(abs(z).'*abs(pr.trg.x))^2 ;
    titer = nan(pr.maxiter,1) ;

    for iter=1:pr.maxiter
        tic
        
        for tt = 1:pr.m
            J(tt, :) = z'*pr.A(tt,:)'*pr.A(tt,:);
        end
        measz = abs(pr.A*z).^2 ;
        F = measz  - pr.y ; 
        Asub_mid = [J, conj(J)];
        %      z_mid = pinv(Asub_mid)*F; % Si Soucis de conditionnement
        %      z_mid = lsqr(Asub_mid,F,1e-2,10);
        beta = 1e-2;
        z_mid = (Asub_mid'*Asub_mid + beta * eye(size(Asub_mid,2)))\(Asub_mid'*F) ;
        z = z - z_mid(1:pr.n); 
        
        titer(iter) = toc ;
        Q(iter+1) = abs(z.'*conj(pr.trg.x))^2/(abs(z).'*abs(pr.trg.x))^2 ;

%         if norm(F)/norm(pr.y) <= pr.tol, break ; end
    end

    zEnd = z ;

    if pr.verbose
        disp([char(9) pr.algo ' : n=' num2str(pr.n) ', m=' num2str(pr.m) ...
                ', epsilon=' num2str(pr.epsilon) ', iter=' num2str(iter) ...
                ', relres=' num2str(norm(abs(measz).^2-pr.y)/norm(pr.y)) ...
                ', relerr=' num2str(norm(pr.trg.x - exp(-1i*angle(trace(pr.trg.x'*z))) * z, 'fro')/norm(pr.trg.x,'fro')) ...
                ', meanIter=' num2str(mean(titer)) 's' ] )
    end 
            
end