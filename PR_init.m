function [ trg, y, A, z0, normest ] = PR_init( pr )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    %% Transfer matrix generation
    switch pr.matrix
        case 'rand'
            A = rand(pr.m,pr.n) ;
        case 'complex_rand'
            A = rand(pr.m,pr.n).*exp(1i*2*pi*rand(pr.m,pr.n)) ;
        case 'custom'
%             load('A16f500_64det46-5px_ideal.mat');
            error('Add custom code here !')
        otherwise
            error('Please choose a correct matrix type : ''complex_rand'', ''rand'', or ''custom'' !')
    end
    A = A/max(max(abs(A))) ; % Normalization
    
    reAbruit = real(A).*(1 + pr.stdBruitMat*randn(size(A))) ;
    imAbruit = imag(A).*(1 + pr.stdBruitMat*randn(size(A))) ;
    Aopt = reAbruit + 1i*imAbruit ;
    
    %% Target values
    trg.x = rand(pr.n,1).*exp(1i*2*pi*rand(pr.n,1)) ; % Target field
    trg.y = (abs(Aopt*trg.x).^2).*(1 + pr.stdBruitMes*randn(size(A,1),1)) ; % Target measurements with noise
    normest = sqrt(sum(trg.y)/numel(trg.y)) ;
    y = trg.y ;
    
    %% Vector to retrieve initialization
    switch pr.initialization
        case 'zeros'
            z0 = ones(pr.n,1) ;
            
        case 'random'
            z0 = 1/sqrt(2)*(randn(pr.n,1)+1i*randn(pr.n,1)) ;
            
        case 'ps'
            z0 = randn(pr.n,1);
            z0 = z0/norm(z0,'fro') ;
            for ii = 1 : pr.psiter
              z0 = A'*(trg.y.* (A*z0)) ;
              z0 = z0/norm(z0,'fro') ;
            end
            z0 = z0*normest ;
            
        case 'null'
            % Chen P. (2017) "Phase Retrieval with One or Two Diffraction
            % Patterns by Alternating Projections with the Null Initialization"
            % doi:10.1007/s00041-017-9536-8
            z0 = randn(pr.n,1);
            z0 = z0/norm(z0,'fro') ;
            for ii = 1 : pr.psiter
              z0 = A'*( ones(size(trg.y)).* (A*z0)) ;
              z0 = z0/norm(z0,'fro') ;
            end
            z0 = z0*normest ;
            
        case 'wfi'
            [nlin,ncol] = size(A);
            lbda = 0;
            for j=1:nlin, lbda = lbda+norm(A(j,:))^2; end
            lbda = sqrt(ncol*sum(trg.y(:))/lbda);
            [lvp,~] = eigs((A'*diag(trg.y(:))*A)/nlin,1);
            z0 = lbda*lvp;
            
        case 'custom'
            error('Add custom code here !')
            
        otherwise
            error('Invalid pr.init parameter, choose among ''zeros'', ''random'', ''ps'', or ''custom'' !')
    end
    
    if ~pr.algoCompare
        zp = reshape(trg.x,sqrt(numel(trg.x))*[1 1]) ;
        figure(1),subplot(2,4,1),imagesc(abs(zp/max(max(abs(zp)))).^2),colorbar,caxis([0 1]),title('Normalized initial intensity')
        figure(1),subplot(2,4,5),imagesc(angle(exp(1i*(angle(zp)-angle(zp(1)))))),colorbar,caxis(pi*[-1 1]),title('Initial phases')
    end
    
end

