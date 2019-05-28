function [ Q, zEnd ] = pr_SSPR( pr )
% Phase retrieval using the SSPR algorithm

    z = pr.z0;
    
    Q(1) = abs(z.'*conj(pr.trg.x))^2/(abs(z).'*abs(pr.trg.x))^2 ;
    titer = nan(pr.maxiter,1) ;    

    u0 = 60 ;
    u = u0 ;
    y1 = 0.9 ;
    y0 = 0.08 ;
    mu = 0.15/1000 ;

    for iter=1:pr.maxiter
        tic
        
        tt = mod(iter,pr.m)+1;
        measz = pr.A*z;
        grad = compute_grad(z,u,sqrt(pr.y),pr.A,tt);
        z = z-mu*grad;
        
        titer(iter) = toc ;
        Q(iter+1) = abs(z.'*conj(pr.trg.x))^2/(abs(z).'*abs(pr.trg.x))^2 ;
        
        if norm(compute_grad(z,u,sqrt(pr.y),pr.A,tt),'fro')< y0*u
            u = y1*u;
        end

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