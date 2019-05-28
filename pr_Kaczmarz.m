function [ Q, zEnd ] = pr_Kaczmarz( pr )
% Phase retrieval using the Kaczmarz algorithm
%
% Bibliography :
%   * K. Wei, "Solving systems of phaseless equations via Kaczmarz methods:
%   a proof of concept study" (2015)

    z = pr.z0;
    
    Q(1) = abs(z.'*conj(pr.trg.x))^2/(abs(z).'*abs(pr.trg.x))^2 ;
    titer = nan(pr.maxiter,1) ;

    for iter=1:pr.maxiter
        maxrelres = 0 ;
        
        tic
        
        for r=1:pr.m
            k=r ;
            k=randi(pr.m) ;
            nrm2 = norm(pr.A(k,:))^2;
            measz = pr.A(k,:)*z;

            maxrelres = max(maxrelres,abs(pr.y(k)-abs(measz)^2)/pr.y(k)) ; 

            z = z + (measz/abs(measz)*sqrt(pr.y(k))-measz)*pr.A(k,:)'/nrm2 ;
        end
        
        titer(iter) = toc ;
        Q(iter+1) = abs(z.'*conj(pr.trg.x))^2/(abs(z).'*abs(pr.trg.x))^2 ;

%         if maxrelres < pr.tol, break ; end
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