%{
    Script combining several phase retrieval methods for comparison

    Jérémy Saucourt, 13/11/2018
%}

clear all
rng('shuffle') % Seeds the random number generator based on current time


%% Problem parameters
%%% Problem size
pr.n = 36 ; % Number of phases to retrieve
pr.m = 8*pr.n ; % Number of detectors
pr.implementedInits = {'zeros','random','ps','null','custom'} ; % 'ps' = power spectrum
pr.initialization = pr.implementedInits{4} ;
pr.matrix = 'complex_rand' ; % Choose between 'complex_rand' (complex), 'rand' (real) and 'custom'
pr.epsilon = 0.0001 ;


%% Phase retrieval algorithm settings
pr.implementedAlgorithms = {'GerchbergSaxton','AlternatedProjections','AlternatedProjectionsPaul', ...
                            'AltMin','SSPR','WirtingerFlow','GaussNewton', ...
                            'Kaczmarz','BlockKaczmarz'} ;
pr.algo = pr.implementedAlgorithms{1} ; % Choose phase retrieval algorithm type
pr.algoCompare = 0 ; % Enable for comparing algorithm performances
pr.maxiter = 20 ; % Maximum internal iterations
pr.psiter = 50 ; % Number of power spectrum iterations (for pr.initialization=='ps')
pr.tol = pr.epsilon ; % Algorithm stops if difference between 2 iterations is <= pr.tol
pr.avg = 10 ; % Averaging
pr.verbose = 1 ; % Enable 'verbose' option to see details
pr.stdBruitMes = 0/100 ; % Ecart-type bruit de mesure d'intensité
pr.stdBruitMat = 0/100 ; % Ecart-type bruit de la matrice de transfert (partie réelle et imaginaire)


%% Phase retrieval
if ~pr.algoCompare % pr.algoCompare == 0 : Single algorithm is tested
    Q = nan(pr.maxiter+1,pr.avg) ;

    for i=1:pr.avg
        tic
        [pr.trg, pr.y, pr.A, pr.z0, pr.normest] = PR_init(pr) ; % Initialization function
        [Q(:,i),zEnd] = PR_algo(pr) ; % Add 'verbose' option to see details
        disp(['moy ' num2str(i) '/' num2str(pr.avg)])
        toc
    end
else % pr.algoCompare == 1 : All algorithms are compared
    Q = nan(pr.maxiter+1,pr.avg,length(pr.implementedAlgorithms)) ;
    
    for i=1:pr.avg
        [pr.trg, pr.y, pr.A, pr.z0, pr.normest] = PR_init(pr) ;
        
        for a=1:length(pr.implementedAlgorithms)
            pr.algo = pr.implementedAlgorithms{a} ;
            [Q(:,i,a),zEnd] = PR_algo(pr) ; % Add 'verbose' option to see details
%             disp([ char(9) 'Avg : ' num2str(i) '/' num2str(pr.avg)])
        end
    end
end


%% Plots
if ~pr.algoCompare % pr.algoCompare == 0 : Single algorithm is tested
    zEnd = reshape(zEnd,sqrt(numel(zEnd))*[1 1]) ;
    figure(1),subplot(2,4,2),imagesc(abs(zEnd/max(max(abs(zEnd)))).^2),colorbar,caxis([0 1]),title('Normalized recovered intensity')
    figure(1),subplot(2,4,6),imagesc(angle(exp(1i*(angle(zEnd)-angle(zEnd(1)))))),colorbar,caxis(pi*[-1 1]),title('Recovered phases')
    
    zInit = reshape(pr.trg.x,sqrt(numel(pr.trg.x))*[1 1]) ;
    figure(1),subplot(2,4,7),imagesc(abs(zEnd).^2-abs(zInit).^2),colorbar,title('Intensity error')
    figure(1),subplot(2,4,8),imagesc(angle(exp(1i*(angle(zEnd)-angle(zEnd(1))-angle(zInit)+angle(zInit(1)))))),colorbar,title('Phase error')

    figure(1)
    subplot(2,4,[3 4]),cla
    plot(Q,'Color','b')
    grid on, box on
    axis([0 pr.maxiter 0 1])
    xlabel('Iteration #')
    ylabel('Recovery phasing quality')

    hold on
    hlim = plot([0 pr.maxiter],[0.99 0.99],'-.g') ;

    if pr.avg > 1
        meanQ = mean(Q') ;
        stdQ = std(Q') ;
        meanStdPQ = meanQ+stdQ  ;
        meanStdMQ = meanQ-stdQ ;
        iter = (1:pr.maxiter)' ;

        hpm = plot(meanQ,'Color','red','LineWidth',3) ;
    end
    % legend([hp(1) hpm hpmp hpmm hlim],{[num2str(algo.nbDisplay) ' tirages'],'Moyenne \mu','\mu + \sigma','\mu - \sigma','\lambda/30 RMS'},'Location','SouthEast')


else % pr.algoCompare == 1 : All algorithms are compared

    linecolors = {[0.9290,0.6940,0.1250],[0.8500, 0.3250, 0.0980],'k',[0,0.5,0],'m',[0.5,0,0.5],'r','c',[0,0,1]} ;
    
    figure(2),cla
        hold on,grid on, box on
        axis([0 pr.maxiter 0 1])
        xlabel('Iteration #')
        ylabel('Recovery phasing quality')
        title('Comparison of several PR algorithms')
    
    if pr.avg == 1
        for a=1:length(pr.implementedAlgorithms)
            h(a) = plot(0:pr.maxiter,Q(:,1,a),'Linewidth',3,'Color',linecolors{a}) ;
        end
        legend([h],pr.implementedAlgorithms,'Location','SouthEast')

    else % pr.avg > 1
        meanQ = squeeze(mean(Q,2)) ;
        stdQ = squeeze(std(Q,0,2)) ;
        meanStdPQ = meanQ+stdQ  ;
        meanStdMQ = meanQ-stdQ ;
        
        for a=1:length(pr.implementedAlgorithms)
            h(a) = plot(0:pr.maxiter,meanQ(:,a),'Linewidth',3,'Color',linecolors{a}) ;
        end
        legend([h],pr.implementedAlgorithms,'Location','SouthEast')
        title(['Average phasing quality (' num2str(pr.avg) ' draws)'])

    end

end

% name = ['Qcomp_n=' num2str(pr.n) '_m,n=' num2str(pr.m/pr.n) '_brMes=' num2str(pr.stdBruitMes*100) '_brMat=' num2str(pr.stdBruitMat*100) '_init=' pr.initialization '_avg=' num2str(pr.avg)] ;
% saveas(figure(2),[pwd '\results\' name '.png'])
