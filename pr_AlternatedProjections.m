function [ Q, zEnd ] = pr_AlternatedProjections( pr )
% Phase retrieval using an Alternated Projections algorithm

    Ainv = pinv(pr.A);
    z = pr.z0;
    
    Q(1) = abs(z.'*conj(pr.trg.x))^2/(abs(z).'*abs(pr.trg.x))^2 ;
    titer = nan(pr.maxiter,1) ;

    for iter=1:pr.maxiter
        tic
        
        measz = pr.A*z;
        yza = sqrt(pr.y).*(measz./abs(measz));
        z = Ainv*(sqrt(pr.y).*exp(1i*angle(pr.A*z))) ;
        
        titer(iter) = toc ;
        Q(iter+1) = abs(z.'*conj(pr.trg.x))^2/(abs(z).'*abs(pr.trg.x))^2 ;

%         if norm(abs(measz).^2-pr.y)/norm(pr.y) <= pr.tol, break ; end
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
