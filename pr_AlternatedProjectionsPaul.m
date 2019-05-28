function [ Q, zEnd ] = pr_AlternatedProjectionsPaul( pr )
% Phase retrieval using an Alternated Projections algorithm

    [U,S,V] = svd(pr.A,0) ;
    s = diag(S) ;

    z = pr.z0 ;
    y = pr.A*z ;
    
    Q(1) = abs(z.'*conj(pr.trg.x))^2/(abs(z).'*abs(pr.trg.x))^2 ;
    titer = nan(pr.maxiter,1) ;
    
    for iter=1:pr.maxiter
        tic
        
        measz = pr.A*z ;
        y = sqrt(pr.y).*exp(1i*angle((U*(U'*y)))) ;
        z = V*((U'*y)./s);
        
        titer(iter) = toc ;
        Q(iter+1) = abs(z.'*conj(pr.trg.x))^2/(abs(z).'*abs(pr.trg.x))^2 ;

%         if pr.algo.critIntStat
%             criter = abs((y'*yanc))^2/(abs(y)'*abs(yanc))^2 ;
%             if criter >= pr.algo.critIntVal, break, end
%         end

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
