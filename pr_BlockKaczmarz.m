function [ Q, zEnd ] = pr_BlockKaczmarz( pr )
% Phase retrieval using the Block-Kaczmarz algorithm
%
% Biliography :
%   * D. Needell, "Randomized Block Kaczmarz Method with Projection for
%   Solving Least Squares" (2014)

    blockSS = round(pr.n/4) ;

    for bls=1:length(blockSS)
        blockSize = blockSS(bls) ;
        nblocks = ceil(pr.m/blockSize) ;
        z = pr.z0 ;

        Q(1) = abs(z.'*conj(pr.trg.x))^2/(abs(z).'*abs(pr.trg.x))^2 ;
        titer = nan(pr.maxiter,1) ;

        for iter=1:pr.maxiter
            maxrelres = 0 ;
            tic

            for r=1:nblocks
                Asub = pr.A((r-1)*blockSize+1:min(r*blockSize,pr.m),:);
                ysub = pr.y((r-1)*blockSize+1:min(r*blockSize,pr.m));
                Asubzbk = Asub*z;

                measz = abs(Asubzbk) ;

                maxrelres = max(maxrelres,norm(ysub-abs(Asubzbk).^2)/norm(ysub));

                z = z +pinv(Asub)*((Asubzbk./abs(Asubzbk)).*sqrt(ysub)-Asubzbk);
            end
            titer(iter) = toc ;

            Q(iter+1) = abs(z.'*conj(pr.trg.x))^2/(abs(z).'*abs(pr.trg.x))^2 ;

%             if maxrelres < pr.tol, break ; end
        end
    end

    zEnd = z ;

    if pr.verbose
        disp([char(9) pr.algo ' : n=' num2str(pr.n) ', m=' num2str(pr.m) ...
            ', epsilon=' num2str(pr.epsilon) ', iter=' num2str(iter) ...
            ', relerr=' num2str(norm(pr.trg.x - exp(-1i*angle(trace(pr.trg.x'*z))) * z, 'fro')/norm(pr.trg.x,'fro')) ...
            ', meanIter=' num2str(mean(titer)) 's' ] )
    end

end