function [ Q, zEnd ] = pr_WirtingerFlow( pr )
% Phase retrieval using the Wirtinger Flow algorithm

    tau0 = 1e8 ; %330
    mu = @(t) min(1-exp(-t/tau0), 1e9) ; %0.2
    z = pr.z0 ;
    
    Q(1) = abs(z.'*conj(pr.trg.x))^2/(abs(z).'*abs(pr.trg.x))^2 ;
    titer = nan(pr.maxiter,1) ;    

    for iter=1:pr.maxiter
        tic
        
        measz = pr.A*z;
        diff = abs(measz).^2-pr.y ;
        grad  = 1/pr.m* pr.A'*( diff .* measz ) ;
        z = z - mu(iter)/pr.normest^2 * grad ;
        
        titer(iter) = toc ;
        Q(iter+1) = abs(z.'*conj(pr.trg.x))^2/(abs(z).'*abs(pr.trg.x))^2 ;
        
        if norm(diff)/norm(pr.y) <= pr.tol, break ; end
    end

    zEnd = z ;

    % Output results of Wirtinger-Flow
    if pr.verbose
        disp([char(9) pr.algo ' : n=' num2str(pr.n) ', m=' num2str(pr.m) ...
                ', epsilon=' num2str(pr.epsilon) ', iter=' num2str(iter) ...
                ', relres=' num2str(norm(abs(measz).^2-pr.y)/norm(pr.y)) ...
                ', relerr=' num2str(norm(pr.trg.x - exp(-1i*angle(trace(pr.trg.x'*z))) * z, 'fro')/norm(pr.trg.x,'fro')) ...
                ', meanIter=' num2str(mean(titer)) 's' ] )
    end
    
end