% [vmatch,H] = MatchGeomVerif(type,f1,f2,match,options) - Geometric verification of image matches
%
% type      = {'PRACT'  ... geometricVerification of Vedaldi's "Practical"
%              'MPRACT'}... geometricVerification + muptiple homographies
% f1, f2    = m x n fature matrix (m=3...disc, 4...oriented disc, 5...ellipse, 6 ... oriented ellipse)
%             as in geometricVerification [[x;y;
% match     = 2 x k matrices [[i1;i2], ...] of indices i1, i2 of matching pairs f1(:,i1), f1(i2,:)
% vmatch    = cell of n 1 x l<=k matrices [j1,i2, ...] of indices j1, ... of verified matches match(:,vmatch{j})
% H         = {H1, H2, ..., Hn} ... homographies of the models
% options   = options {'name1',value1,'name2','value2',...}

% 2015-05-14 Tomas Pajdla, pajdla@cmp.felk.cvut.cz
%            adapted from A. Vedaldi's geometricVerification.m in "Practical"   
function [Ins,H21] = MatchGeomVerif(type,f1,f2,matches,varargin)

switch type
    case 'PRACT'
        % geometricVErification.m default values
        opts.tolerance1 = 20;
        opts.tolerance2 = 20;
        opts.tolerance3 = 10;
        opts.minInliers = 6;
        opts.numRefinementIterations = 3;
        opts.stopInlierFraction = 0.7;
        opts = vl_argparse(opts, varargin) ;
        
        numMatches = size(matches,2) ;
        inliers = cell(1, numMatches) ;
        
        x1 = double(f1(1:2, matches(1,:))) ;
        x2 = double(f2(1:2, matches(2,:))) ;
        
        x1hom = x1 ;
        x2hom = x2 ;
        x1hom(end+1,:) = 1 ;
        x2hom(end+1,:) = 1 ;
        
        % bad set of candidate inliers will produce a bad model, but
        % this will be discared
        warning('off','MATLAB:rankDeficientMatrix') ;
        
        for m = 1:numMatches
            for t = 1:opts.numRefinementIterations
                if t == 1
                    A1 = toAffinity(f1(:,matches(1,m))) ;
                    A2 = toAffinity(f2(:,matches(2,m))) ;
                    H21 = A2 * inv(A1) ;
                    x1p = H21(1:2,:) * x1hom ;
                    tol = opts.tolerance1 * sqrt(det(H21(1:2,1:2))) ; 
                elseif t <= 4
                    % affinity
                    H21 = x2(:,inliers{m}) / x1hom(:,inliers{m}) ;
                    x1p = H21(1:2,:) * x1hom ;
                    tol = opts.tolerance2 * sqrt(det(H21(1:2,1:2))) ;
                else
                    % homography
                    x1in = x1hom(:,inliers{m}) ;
                    x2in = x2hom(:,inliers{m}) ;
                    
                    % Sanity check
                    %H = [.1 0 .4 ; 2 .3 .5 ; .1 .002 1] ;
                    %x1in = [randn(2,100) ; ones(1,100)] ;
                    %x2in = H*x1in ;
                    %x2in = bsxfun(@times, x2in, 1./x2in(3,:)) ;
                    
                    S1 = centering(x1in) ;
                    S2 = centering(x2in) ;
                    x1c = S1 * x1in ;
                    x2c = S2 * x2in ;
                    
                    M = [x1c, zeros(size(x1c)) ;
                        zeros(size(x1c)), x1c ;
                        bsxfun(@times, x1c,  -x2c(1,:)), bsxfun(@times, x1c,  -x2c(2,:))] ;
                    [H21,D] = svd(M,'econ') ;
                    H21 = reshape(H21(:,end),3,3)' ;
                    H21 = inv(S2) * H21 * S1 ;
                    H21 = H21 ./ H21(end) ;
                    
                    x1phom = H21 * x1hom ;
                    x1p = [x1phom(1,:) ./ x1phom(3,:) ; x1phom(2,:) ./ x1phom(3,:)] ;
                    tol = opts.tolerance3 * sqrt(det(H21(1:2,1:2))) ;
                end
                dist2 = sum((x2 - x1p).^2,1) ;
                inliers{m} = find(dist2 < tol^2) ;
                if numel(inliers{m}) < opts.minInliers, break ; end
                if numel(inliers{m}) > opts.stopInlierFraction * size(matches,2), break ; end % enoguh!
            end
        end
        scores = cellfun(@numel, inliers) ;
        [~, best] = max(scores) ;
        Ins = {inliers{best}} ;
    case 'MPRACT' % Practical's geometric verification with multiple homographies
        % input parameters f1,f2,matches,varargin
        % set default options
        opts.SimTol = 10;
        opts.AffTol = 10;
        opts.HomTol = 5;
        opts.MinIns = 6;
        opts.RefIterNum = 8;
        opts.StopInsFrac = 0.7;
        opts.MaxHoms = 3;
        opts.MinInsNum = 20;
        % inint output
        Ins = {[]};
        % reset specified options
        opts = optsparse(opts, varargin);
        % init matches
        ms = matches;
        mi = 1:size(matches,2); % match index
        for i=1:opts.MaxHoms % up to maximal number of homograpjies per image pair
            % init inliers
            ins = cell(1,size(ms,2));
            % coordinates of matches
            x1 = a2h(double(f1(1:2,ms(1,:))));
            x2 = a2h(double(f2(1:2,ms(2,:))));
            % bad set of candidate inliers will produce a bad model, but this will be discared
            warning('off','MATLAB:rankDeficientMatrix') ;
            % initialize and grow from all input matches
            for m = 1:size(ms,2)
                for t = 1:opts.RefIterNum
                    if t == 1
                        % similarity
                        A1 = toAffinity(f1(:,ms(1,m)));
                        A2 = toAffinity(f2(:,ms(2,m)));
                        H21 = A2*inv(A1);
                        x1p = H21(1:2,:)*x1;
                        tol = opts.SimTol; % *sqrt(det(H21(1:2,1:2))); % may be unstable
                    elseif t <= 4
                        % affinity
                        H21 = x2(1:2,ins{m})/x1(:,ins{m});
                        x1p = H21(1:2,:)*x1;
                        tol = opts.AffTol; % *sqrt(det(H21(1:2,1:2))); % may be unstable
                    else
                        % homography
                        x1in = x1(:,ins{m});
                        x2in = x2(:,ins{m});
                        S1 = centering(x1in);
                        S2 = centering(x2in);
                        x1c = S1*x1in;
                        x2c = S2*x2in;
                        M = [                         x1c               zeros(size(x1c))
                            zeros(size(x1c))                            x1c
                            bsxfun(@times,x1c,-x2c(1,:))  bsxfun(@times,x1c, -x2c(2,:))];
                        [H21,~] = svd(M,'econ');
                        H21 = reshape(H21(:,end),3,3)';
                        H21 = inv(S2)*H21*S1;
                        H21 = H21/H21(end);
                        x1phom = H21*x1;
                        x1p = h2a(x1phom); % [x1phom(1,:)./x1phom(3,:) ; x1phom(2,:) ./ x1phom(3,:)] ;
                        tol = opts.HomTol; % *sqrt(det(H21(1:2,1:2))); % unstable
                    end
                    % core the model
                    d2 = sum((x2(1:2,:)-x1p).^2,1);
                    ins{m} = find(d2<tol^2);
                    if m==146
                        m;
                    end
                    if numel(ins{m})<opts.MinIns, break ; end
                    if numel(ins{m})>opts.StopInsFrac*size(ms,2), break ; end % enoguh!
                end
            end
            % select the best model
            scores = cellfun(@numel,ins);
            [~, best] = max(scores);
            bin = ins{best}; % inliers of the best model
            if numel(bin)<opts.MinInsNum, break; end % STOP when the models get too small                
            Ins{i} = mi(bin); % store the index of the verified matches, indexing the input list of matches 
            ix = setdiff(1:size(ms,2),bin); % remaining tentative maches index
            if numel(ix)<opts.MinInsNum, break; end % STOP when the number of remaining matches is too small
            ms = ms(:,ix); % remaining tentative matches
            mi = mi(ix); % index of the remaining tentative natches into the input list of matches
        end
    otherwise
        error('MatchGeomVerif: %s type not implemented.',type);
end
end

% --------------------------------------------------------------------
function C = centering(x)
% --------------------------------------------------------------------d
T = [[1 0;0 1] -mean(x(1:2,:),2); 0 0 1];
sx1 = std(x(1,:)); if sx1<0.1, sx1=0.1; end
sx2 = std(x(2,:)); if sx2<0.1, sx2=0.1; end
x = T*x;
S = [1/sx1 0     0
         0 1/sx2 0 
         0     0 1];
C = S*T;
end

% --------------------------------------------------------------------
function A = toAffinity(f)
% --------------------------------------------------------------------
switch size(f,1)
    case 3 % discs
        T = f(1:2) ;
        s = f(3) ;
        th = 0 ;
        A = [s*[cos(th) -sin(th) ; sin(th) cos(th)], T ; 0 0 1] ;
    case 4 % oriented discs
        T = f(1:2) ;
        s = f(3) ;
        th = f(4) ;
        A = [s*[cos(th) -sin(th) ; sin(th) cos(th)], T ; 0 0 1] ;
    case 5 % ellipses
        T = f(1:2) ;
        A = [mapFromS(f(3:5)), T ; 0 0 1] ;
    case 6 % oriented ellipses
        T = f(1:2) ;
        A = [f(3:4), f(5:6), T ; 0 0 1] ;
    otherwise
        assert(false) ;
end
end

% --------------------------------------------------------------------
function A = mapFromS(S)
% --------------------------------------------------------------------
% Returns the (stacking of the) 2x2 matrix A that maps the unit circle
% into the ellipses satisfying the equation x' inv(S) x = 1. Here S
% is a stacked covariance matrix, with elements S11, S12 and S22.

tmp = sqrt(S(3,:)) + eps ;
A(1,1) = sqrt(S(1,:).*S(3,:) - S(2,:).^2) ./ tmp ;
A(2,1) = zeros(1,length(tmp));
A(1,2) = S(2,:) ./ tmp ;
A(2,2) = tmp ;
end