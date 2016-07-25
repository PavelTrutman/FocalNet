% [p,inl,res,L] = ransacfit(x,model,maxres,smpls,smpln[,doLO,cnstr]) - Ransac fitting
%
% x         ... data - structure depends on the model, see the function
%                      implementation
% model     ... {'line','line-cnstr-seg','plane','conic','3Dpts','E5ptNister','P3PSolver','uX2PSolver'}
% maxres    ... maximal residual to vote
% smpls     ... sample size 0    = minimal sample
%                          (0,1) = percentile from the size of data
%                           >=1  = number of point in the samples
% smpln     ... number of samples
% doLO      ... {0,1} do local optimization (0 implicit)
% cnstr     ... constraint parameters
%               'line-cnstr-seg' ... cnstr = [x1 x2], such that x1 resp. x2 
%                                                     is on the left resp. right 
%                                                     from the line 
% p         ... model parameters
% inl       ... inliers
% res       ... residuals for the model p
% L         ... fitting structure L{i} = output at the i-th improvement
%
% Example: >>ransacfit;

% T.Pajdla, pajdla@cmp.felk.cvut.cz
% 2015-07-11
function [p,inl,r,L] = ransacfit(x,model,maxres,smpls,smpln,doLO,cnstr)
 if nargin>0 %% Fit
     if nargin<7, cnstr = []; end
     if nargin<6, doLO = false; end % implicitly do not use the local optimization
     %% model independent sample size
     % percentile from the data lengths
     if smpls >0 && smpls <1, smpl = max(2,ceil(length(x)*smpls)); end
    switch model
        case 'line'
            if smpls == 0, smpl = 2; end % the minimal case
            if smpls>=2, smpl = smpls; end  % enforce the number of samples
            fitF = @l2Dfit; % line fit function
            resF = @l2Dres; % residual function 
            fitDataF = @ransacfitIdF; % original data
            resDataF = @ransacfitIdF; % original data            
        case 'line-cnstr-seg'
            if smpls == 0, smpl = 2; end % the minimal case
            if smpls>=2, smpl = smpls; end  % enforce the number of samples
            fitF = @l2Dfit; % line fit function
            resF = @l2DresCnstrSeg; % residual function with constraints enforced
            fitDataF = @ransacfitIdF; % original data
            resDataF = @ransacfitIdF; % original data    
        case 'plane'
            if smpls == 0, smpl = 3; end % the minimal case
            if smpls>=3, smpl = smpls; end  % enforce the number of samples
            fitF = @p3Dfit; % line fit function
            resF = @p3Dres; % residual function 
            fitDataF = @ransacfitIdF; % original data
            resDataF = @ransacfitIdF; % original data                        
        case 'conic'
            if smpls == 0, smpl = 5; end % the minimal case
            if smpls>=5, smpl = smpls; end  % enforce the number of samples
            fitF = @fitellipse;
            resF = @elliRes; % residual function   
            fitDataF = @ransacfitIdF; % original data
            resDataF = @ransacfitIdF; % original data           
        case '3Dpts'
            if smpls == 0, smpl = 3; end % the minimal case
            if smpls>=3, smpl = smpls; end  % enforce the number of samples
            fitF = @XY2rts; % best similarity fit of two sets of 3D points
            resF = @XY2rtsErr; % residual function  
            fitDataF = @ransacfitIdF; % original data
            resDataF = @ransacfitIdF; % original data
        case 'E5ptNister'
            % x is struct containing
            % x.iK = [iK1;iK2] ... stacked inversions of camera calibration matrices
            % x.x = [x1;x2] ... stacked homogenenous coordinates of tentative image matches
            if smpls == 0, smpl = 5; end % the minimal case
            if smpls>=5, smpl = smpls; end  % enforce the number of samples
            fitF = @E5ptNister; % Nister's 5pt method
            resF = @EGeomErr; % residual function
            fitDataF = @E5ptFitImPoints; % compute calibrated image measurements
            resDataF = @E5ptResImPoints; % pass the original uncalibrated data
            cnstr = x.iK; % store inverses of the camera calibration matrices for computing F in EGeomErr
        case 'P3PSolver'
            % x is struct containing
            % x.K = camera calibration matrix
            % x.iK = inversion of the camera calibration matris
            % x.x = [u;X] ... stacked hom coordis of image points & coords of 3D points
            if smpls == 0, smpl = 3; end % the minimal case
            if smpls>=3, smpl = smpls; end  % enforce the number of samples
            fitF = @P3PSolver; % 4pt absolute pose
            resF = @PerspRepErr; % residual function            
            fitDataF = @P3PSolverFitImPoints; % original data
            resDataF = @E5ptResImPoints; % original data   
            cnstr = x.K; % store camera calibration matrix for computing P in PerspRepErr
        case 'uX2PSolver'
            % x is struct containing
            % x.K = camera calibration matrix
            % x.iK = inversion of the camera calibration matris
            % x.x = [u;X] ... stacked hom coordis of image points & coords of 3D points
            if smpls == 0, smpl = 6; end % the minimal case
            if smpls>=6, smpl = smpls; end  % enforce the number of samples
            fitF = @uX2PSolver; % DLT absolute pose
            resF = @PerspRepErr; % residual function            
            fitDataF = @P3PSolverFitImPoints; % original data
            resDataF = @E5ptResImPoints; % original data   
            cnstr = x.K; % store camera calibration matrix for computing P in PerspRepErr
        otherwise
            error('ransacfit: unknown model');
    end
    % check if the sample size has been assigned
    if ~exist('smpl','var'), error('ransacfit: invalid sample size for smpls = %s',smpls); end            
    %% RANSAC
    % prepare method depndent data from fitting and evaluation
    xf = fitDataF(x);
    xr = resDataF(x);
    % do ransac
    si = round((size(xf,2)-1)*rand(smpl,smpln)+1); % samples in the columns
    iN = 0;
    p = [];
    inl = [];
    r = [];
    k = 1; 
    for i = si
        y = xf(:,i); % get a sample
        pf = fitF(y); % fit a model
        r  = resF(pf,xr,cnstr); % eval residuals & constraints
        il = abs(r)<=maxres; % inliers
        if size(il,1)>1 % there are more alternative models returned by fitF
            in = sum(il,2); % inlier # 
            [in,ix] = max(in);
            il = il(ix,:); % select the best inliers
            r = r(ix,:); % select the residuals
            pf = pf{ix}; % select the best model
        else
            in = sum(il); % inlier #
            if ~isempty(pf)
                if iscell(pf)
                    pf = pf{1};
                end
            end
        end
        if in>iN % larger support
            if doLO % local optimization
                y = xf(:,il); % get all inliers
                po = fitF(y); % fit a model
                ro = resF(po,xr,cnstr); % eval residuals & constraints
                ilo = abs(ro)<=maxres; % inliers
                if size(ilo,1)>1 % there are more alternative models returned by fitF
                    ino = sum(ilo,2); % inlier #
                    [ino,ix] = max(ino);
                    ilo = ilo(ix,:); % select the best inliers
                    ro = ro(ix,:); % select the residuals
                    po = po{ix}; % select the best model
                else
                    ino = sum(ilo); % inlier #
                end
                if ino>in % an improvement
                    p = po; r = ro; il = ilo; in = ino; 
                end
            end
            iN  = in;
            inl = il;
            p = pf;
            if nargout>3 % store the history
                L{k}.iN = iN; L{k}.p = p; L{k}.res = r; L{k}.inl = inl; 
                k = k + 1;
            end
        end
    end
    % refit to all inliers (must return only one model!)
    if sum(inl)>0
        y   = xf(:,inl); % inliers
        p = fitF(y,p); % fit
        r = resF(p,xr,cnstr); % residuals
        inl = abs(r)<=maxres; % inliers
    else
        
    end
    if nargout>3 % store the history
        L{k}.iN = iN; L{k}.p = p; L{k}.res = r; L{k}.inl = inl;
    end
 else  %% Demo
    %randn('state',1); rand('state',1);
    t = 0:100;
    x = [t;1/2*t+1];
    x = x+randn(size(x));
    ixi = [1:20 30:40 70:90];
    x(:,ixi) = x(:,ixi)+10*randn(2,52);
    mr = 1.65;
%    [p,inl,res] = ransacfit(x,'line-cnstr-seg',mr,2,20,0,[0 15;0 20]');
    [p,inl,res] = ransacfit(x,'line',mr,2,20);    
    subfig(2,2,1);
    plot(x(1,inl),x(2,inl),'.','markersize',10); hold
    plot(x(1,~inl),x(2,~inl),'.r','markersize',10)
    axis([-10 110 -10 70]);
    if ~isempty(p),
        set(plotline(p),'linewidth',2);
        title('RANSAC linefit: b - inliers, r - outliers')
        subfig(2,2,2);
        ix = 1:size(x,2);
        plot(ix(inl),res(inl),'.b',ix(~inl),res(~inl),'.r'); hold
        axis([-10 110 -25 25]);
        set(plotline([[0;1;mr],[0;1;-mr]]),'color','g');
        plotline([0;1;0]);
        title('RANSAC residuals: b - inliers, r - outliers')
        subfig(2,2,3);
        plot(x(1,inl),x(2,inl),'.','markersize',10); hold
        axis([-10 110 -10 70]);
        set(plotline(p),'linewidth',1);
        ixn = true(size(ixi));
        ixn(ixi) = false;
        ix = 1:length(ix);
        ix = ix(ixn);
        for i=1:10
            p = x(:,randsample(ix,2));
            if i==2
                l = cross([p(1,1);p(2,1);1],[p(1,2);p(2,2);1]);
                set(plotline(l),'color','g','linewidth',2')
                plot(p(1,:),p(2,:),'w.','markersize',30);
                plot(p(1,:),p(2,:),'g.','markersize',20);
            end
        end
    end
end
    

