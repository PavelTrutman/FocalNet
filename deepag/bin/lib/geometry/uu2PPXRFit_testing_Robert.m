% [P,X,in,e,F,u] = uu2PPXRFit(u,K,op,debug) - calibrated reconstruction from two images
%
% u = {u1, u2} with u1, u2 2 x n image projection matrices
% K = {K1, K2} with 3 x 3 camwra calibration matrices
% op = parameters
%
%      RANSAC (see ransacfit with E5ptNister model):
%      op.maxres = 2*sN;  % maximal residual in pixels = 2*max image noise
%      op.smpls  = 0;     % minimal sample size
%      op.smpln  = 1000;  % number of samples
%      op.doLO   = false; % do Local Optimization
%      op.cnstr  = [];    % constraints
%      op.V2TM   = to wiew tringulation method, see uu2X.m
% P = {P1, P2} with 3 x 4 camera projection matrices,
%              P1 = [I|0] and P2 is selected to maximize the number of
%                         inliers reconstructed infront of cameras
% X = 3 x m reconstructed points for inliers infront of the cameas
% in = inlier reconstructed infront of cameras index into u
% F = fundamental matrix
% e = image reprojection errors for all points
% debug = show debugging plots if true
%
% T. Pajdla, pajdla@cmp.felk.cvut.cz, 2015-09-08
function [P,X,in,e,F,x,f] = uu2PPXRFit(x,K,op,debug)
if nargin>0
    if nargin<4
        debug = false;
    end
    f = []; % no figures will be plot
    % preapare data for E5ptRansacing
    xiK.x  = [a2h(x{1}); a2h(x{2})]; % uncalibrated image points
    xiK.iK = [inv(K{1}); inv(K{2})]; % inverses of camera calibration matrices
    [E,in] = ransacfit(xiK,'E5ptNister',op.maxres,op.smpls,op.smpln,op.doLO,op.cnstr); % run ransac
    if ~isempty(in)
        F = xiK.iK(4:6,:)'*E*xiK.iK(1:3,:); % get F for original image coordinates
        e = EGeomErr(F,xiK.x); % epipolar errors
        P = E2PP(E); % camera matrices with P1 = [I 0]
        P{1} = cellfun(@(x) xiK.iK(1:3,:)\x,P{1},'UniformOutput',false);
        P{2} = cellfun(@(x) xiK.iK(4:6,:)\x,P{2},'UniformOutput',false);
        % select P2 to maximize the number of points reconstruted in front of the cameras
        X = cell(numel(P{1}),1);
        zp = cell(numel(P{1}),1);
        nzp = zeros(numel(P{1}),1);
        for i=1:numel(P{1})
            [X{i},~,z] = uu2X([x{1}(:,in);x{2}(:,in)],[P{1}{i};P{2}{i}],op.V2TM); % 3D points
            zp{i} = all(z>0); % point in front of both cameras and small reprojection error
            nzp(i) = sum(zp{i}); % the number of points infront of both cameras
        end
        [~,mix] = max(nzp); % index of the projection matrices with a maximal number of points
        if debug
            close all;
            for i=1:numel(X)
                subfig(2,4,i); plot3d(X{i},'.'); hold;
                text(X{i}(1,[1 end]),X{i}(2,[1 end]),X{i}(3,[1 end]),num2str(([1 size(X{i},2)])'));
                camplot(P{1}{i}/vnorm(P{1}{i}(3,1:3))/5);[~,~,C1]=P2KRC(P{1}{i}); text(C1(1),C1(2),C1(3),'1','fontsize',15,'color','r');
                camplot(P{2}{i}/vnorm(P{2}{i}(3,1:3))/5);[~,~,C2]=P2KRC(P{2}{i}); text(C2(1),C2(2),C2(3),'2','fontsize',15,'color','r');
                axis equal; grid; title(sprintf('# z>0 = %d',nzp(i)));
            end
        end
        P{1} = P{1}{mix}; % select cameras
        P{2} = P{2}{mix};
        X = X{mix}; % select points
        zp = zp{mix}; % the positive ones
        in(in) = in(in) & zp; % update inliers to be infront of cameras
        X = X(:,zp); % select only inlier 3D points
    else
        P = {nan(3,4),nan(3,4)}; X = []; e = []; F =[]; x =[];
    end
else % unit tests
    % orientation
    X = [ 0 10 10  0 0  10 10  0
          0  0 10 10 0   0 10 10 
         10 10 10 10 15 15 15 15];
    X = bsxfun(@minus,X,[5;5;-5]);
    K = [1000     0   500
              0  1000 500
              0    0   1];
    C1 = [0;0;0];
    C2 = [2;0;0];
     A = [ 1  0 0
           0  1 0
           0  0 1];
    P1 = K*[eye(3) -C1];
    P2 = K*[eye(3) -C2];
    u1 = h2a(A*X2u(X,P1)); 
    u2 = h2a(A*X2u(X,P2)); 
    subfig(2,4,1); plot3d(X,'.'); hold; camplot(P1); camplot(P2); axis equal; grid
    text(X(1,[1:end]),X(2,[1:end]),X(3,[1:end]),num2str([1:size(X,2)]'));
    text(C1(1),C1(2),C1(3),'1','color','b','fontsize',15);
    text(C2(1),C2(2),C2(3),'2','color','b','fontsize',15);
    subfig(2,4,2); plot3d(u1(1:2,:),'.');axis([0 1000 0 1000]);title('cam 1');
    subfig(2,4,3); plot3d(u2(1:2,:),'.');axis([0 1000 0 1000]);title('cam 2');
    op.maxres = 1; % maximal residual in pixels
    op.smpls = 0; % minimal sample size
    op.smpln = 1000; % number of samples
    op.doLO = false; % do Local Optimization
    op.cnstr = []; % constraints 
    op.V2TM = 'TRAN'; % 2D triangulation at 1/2 of the shortest transversal
    [PP,Y,in,e,F] = uu2PPXRFit({u1 u2},{K K},op,true);
    pause
    % simulate a scene
    clear K R C P
    x = trandn([1 500],3,-2,2);
    y = trandn([1 500],3,-2,2);
    z = trandn([1,500],0.5,-0.5,0.5)+7;
    X = [x(:)';y(:)';z(:)'];
    clear x y z
    % simulate cameras with pixel size 0.005 mm and f = 1000 pixels = 5 mm
    K{1} = [1/0.005       0  500/(1000*0.005)
             0      1/0.005  500/(1000*0.005)
             0            0    1/(1000*0.005)];
    K{2} = K{1};
    R{1} = a2r([1;1;1],round(3*randn(1))/10);
    R{2} = a2r([1;1;1],round(3*randn(1))/10);
    C{1} = [-1;0;2];
    C{2} = [ 1;0;-5];
    P{1} = K{1}*[R{1} -R{1}*C{1}];
    P{2} = K{2}*[R{2} -R{2}*C{2}];
    Ft = PP2F(P{1},P{2}); % F contructed from projection matrices
    C2t = P{2}(:,end); % second ground truth camera center
    % projections
    x{1} = h2a(P{1}*a2h(X));
    x{2} = h2a(P{2}*a2h(X));
    % add noise
    sN = 0.5;
    x{1} = x{1}+sN*trandn(size(x{1}),sN,-3*sN,3*sN);
    x{2} = x{2}+sN*trandn(size(x{2}),sN,-3*sN,3*sN);
    % generate mismatches
    N = floor(2*size(X,2)/3);
    % rng('default'); % repeatable random sequences
    oix = randperm(size(X,2),N); % random positions to swap
    oix = [oix; sort(oix)]; % swap index
    oix = oix(:,diff(oix)~=0); % remove fixed points
    x{2}(:,oix(2,:)) = x{2}(:,oix(1,:)); % swap
    % compute the epipolar geometry and select matches
    % preapare data for E5ptRansacing
    op.maxres = 2*sN; % maximal residual in pixels
    op.smpls = 0; % minimal sample size
    op.smpln = 1000; % number of samples
    op.doLO = false; % do Local Optimization
    op.cnstr = []; % constraints 
    op.V2TM = 'TRAN'; % 2D triangulation at 1/2 of the shortest transversal
    % compute the reconstruction
    [PP,Y,in,e,F] = uu2PPXRFit(x,K,op,true);
    pause
    % Normalize PP{2} translation to get the same scale as in the GT
    nt = vnorm(Ft(:))/sqrt(2); % length of the GT camera translation
    PP{2}(:,end) = nt*PP{2}(:,end);
    H = [R{1} -R{1}*C{1};0 0 0 1];
    PP = cellfun(@(x) x*H,PP,'UniformOutput',false);
    % Y = h2a(inv([R{1} -R{1}*C{1};0 0 0 1])*a2h(nt*Y));
    Y = h2a(H\a2h(nt*Y));
    % plot the results  
    f(1) = subfig(2,4,1);
    plot3d(X,'.');hold
    h = cellfun(@plotCam,P,'UniformOutput',false);axis equal; grid;
    cellfun(@(x) set(x,'color','b'),h);
    title('3D points & cameras');
    f(2) = subfig(2,4,2);
    ix = 1:length(e);
    ouTP = intersect(oix(2,:),find(~in)); % True Positive Outliers
    ouFP = setdiff(find(~in),oix(2,:)); % False Positive Outliers
    ouFN = setdiff(oix(2,:),find(~in)); % False Negative Outliers
    ouTN = intersect(setdiff(ix,oix(2,:)),find(in)); % True Negative Outliers
    em = max(e,0.1);
    %semilogy(ix,em,'-r');
    display(oix)
    %display(em(oix(2,:)))
    display('test')
    pause
    semilogy(oix(2,:),em(oix(2,:)),'.r');hold;
    semilogy(setdiff(ix,oix(2,:)),em(setdiff(ix,oix(2,:))),'.g');
    semilogy(ix(~in),em(~in),'or'); 
    semilogy(ix(in),em(in),'ob');
    semilogy(ix([1 end]),[op.maxres op.maxres],'k');
    axis tight; title('F estim: Err = max(d(x2,F*x1),d(x2''*F,x1))');
    f(3) = subfig(2,2,2);
    for i=1:numel(P)
        subplot(1,numel(P),i);
        plot3d(x{i},'.b'); hold; plot3d(x{i}(:,oix(1,:)),'.r'); plot3d(x{i}(:,in),'og'); 
        axis ij; axis equal; axis([0 1000 0 1000]);
        title(sprintf('Camera %d',i));
    end
    f(4) = subfig(2,6,7);
    plottab([numel(ouTP) numel(ouFN);numel(ouFP) numel(ouTN)]/numel(ix),...
        ['TNM';'TPM'],['DNM';'DPM'],'True\\Dcsn','%.3f');
    title('Positive, Negative, Matches');
    f(5) = subfig(2,6,8);
    plottab([Ft;F],['  ';'TF';'  ';'  ';' F';'  '],'','','%.3f');
    title('F mat: TF - true F, F - estim')
    f(6) = subfig(2,6,9);
    TC2 = P{2}(:,1:3)\P{2}(:,end); 
    CC2 = PP{2}(:,1:3)\PP{2}(:,end);
    plottab([TC2',NaN;
             CC2',NaN;
             (CC2-TC2)',vnorm(CC2-TC2)],...
             ['C2t    ';'C2     ';'C2-C2t '],[' ';'v';' ';'d'],'','%.2f');
    title('2nd cam C: C2t - true, C2 - estim');
    figure(f(1))
    plot3d(Y,'.r');
    h = cellfun(@plotCam,PP,'UniformOutput',false);
    cellfun(@(x) set(x,'color','r'),h);    
    axis equal; grid on;
    title('3D points & cams: b=true, r=estim')
    P = PP;
    X = Y;
end

