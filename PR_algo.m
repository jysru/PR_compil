function [ Q, zEnd ] = PR_algo( pr )
% Selects among implemented algorithms for Phase Retireval
%
% Most of the PR algorithms listed here have been extracted from the
% PhasePack package available at:
% http://cs.umd.edu/~tomg/projects/phasepack/
% 
% For additional information, check out the user guide available at:
% https://arxiv.org/abs/1711.09777

    Q = nan(pr.maxiter,1) ;

    switch pr.algo
        case 'GerchbergSaxton'
            [ Q, zEnd ] = pr_GerchbergSaxton( pr ) ;
        case 'AlternatedProjections'
            [ Q, zEnd ] = pr_AlternatedProjections( pr ) ;
        case 'AlternatedProjectionsPaul'
            [ Q, zEnd ] = pr_AlternatedProjectionsPaul( pr ) ;
        case 'SSPR'
            [ Q, zEnd ] = pr_SSPR( pr ) ;
        case 'WirtingerFlow'
            [ Q, zEnd ] = pr_WirtingerFlow( pr ) ;
        case 'AltMin'
            [ Q, zEnd ] = pr_AltMin( pr ) ;
        case 'GaussNewton'
            [ Q, zEnd ] = pr_GaussNewton( pr ) ;
        case 'Kaczmarz'
            [ Q, zEnd ] = pr_Kaczmarz( pr ) ;
        case 'BlockKaczmarz'
            [ Q, zEnd ] = pr_BlockKaczmarz( pr ) ;
        otherwise
            error('Please choose between implemented algorithms : ''gerchberg-saxton'', ''wirtinger-flow'', ''kaczmarz'' or ''block-kaczmarz'' !')
    end

end

