% [H,e] = H4xP2P(P1,P2[,op]) - Matrix H fitting external orientation of cameras P1 to cameras P2 
%
% P1, P2 = celarray (or matrices) of 3x4 rank 3 projection matrix
% H = 4x4 matrix sich that P2{:} ~ P1{:}*H
% e = residual error
% op = options
%      '1stto1st'    = exact fit for the first cameras
%      'fitCfitR' = fit C and R by creating a point on the optical axis of each camera
%                   in the distace of max dist between the camera centers (implicit)
%      'setCfitR' = set C and if any freedom left, fit R
%      'setR&|b|' = set R and the distance between the centers
%      'fitC2C'   = fit camera centers to camera centers

% T.Pajdla, pajdla@cmp.felk.cvut.cz
% 2015-08-26
function [H,e] = H4xP2P(P1,P2,op)
debug = false; %true;
if nargin>0
    if nargin<3, op = 'fitCfitR'; end
    % extract first matrices
    if ~iscell(P1), P1 = {P1}; end;
    if ~iscell(P2), P2 = {P2}; end;
    % decompose P's
    [~,R1,C1] = cellfun(@(x) P2KRC(x), P1, 'UniformOutput', false);
    [~,R2,C2] = cellfun(@(x) P2KRC(x), P2, 'UniformOutput', false);
    switch op
        case '1stto1st' %  first cameras
            % find the transform
            E1 = [R1{1} -R1{1}*C1{1};[0 0 0 1]];
            E2 = [R2{1} -R2{1}*C2{1};[0 0 0 1]];
            H  = E1\E2;
            e = 0;
        case 'setCfitR'
            error('setCfitR not implemented.')
        case 'setR&|b|'
            d1  = vnorm(C1{2}-C1{1}); % baseline length of the first camera pair
            d2  = vnorm(C2{2}-C2{1}); % baseline length of the second camera pair
            E1 = [R1{1} -R1{1}*C1{1};[0 0 0 1]];
            E2 = [R2{1} -R2{1}*C2{1};[0 0 0 1]];
            H  = E2\E1; % point transformation
            H  = [eye(3) -C2{1};0 0 0 1]*H; % C2{1} to origin
            H  = diag([d2/d1*[1 1 1] 1])*H; % scale
            H  = [eye(3)  C2{1};0 0 0 1]*H; % back to C2{1}
            H  = H\eye(4); % camera transformation
            e  = [NaN NaN];
        case 'fitC2C' % fit camera centers to camera centers
            C1 = cat(2,C1{:});
            C2 = cat(2,C2{:});
            H = XY2rts([C1;C2]); % optimal registration
            H = [H.s*H.r H.t;0 0 0 1]; % mapping of points
            e = vnorm(C2-h2a(H*a2h(C1)));
            H = H\eye(4); % mapping of cameras
        case 'fitCfitR'
            c1 = cat(2,C1{:});
            c2 = cat(2,C2{:});
            d1  = mean(vnorm(kron(c1,ones(1,size(c1,2)))-kron(ones(1,size(c1,2)),c1))); % size of the scene 1
            d2  = mean(vnorm(kron(c2,ones(1,size(c2,2)))-kron(ones(1,size(c2,2)),c2))); % size of the scene 1
            V1 = cellfun(@(r,c) r'*[0;0;d1]+c,R1,C1,'UniformOutput',false); % a reasonble point infront of the first set of cameras
            V2 = cellfun(@(r,c) r'*[0;0;d2]+c,R2,C2,'UniformOutput',false); % a reasonble point infront of the second set of cameras
            X1 = [c1 cat(2,V1{:})]; % two sets of points to register
            X2 = [c2 cat(2,V2{:})];
            H = XY2rts([X1;X2]); % optimal registration
            H = [H.s*H.r H.t;0 0 0 1]; % mapping of points
            if debug
                subfig(2,3,1);
                plot3d(X2(:,1:size(X2,2)/2),'ob'); hold; plot3d(X2(:,size(X2,2)/2+1:end),'.b');
                X1H = h2a(H*a2h(X1));
                plot3d(X1H(:,1:size(X2,2)/2),'or'); plot3d(X1H(:,size(X2,2)/2+1:end),'.r');
                grid
            end
            H = H\eye(4); % mapping of cameras
            e  = [NaN NaN];
    end
else % unit tests
    % one camera
    P1 = [1000    0 500 0
             0 1000 500 0
             0    0   1 0];
    Ht  = [a2r([1;1;1],pi/5)*[eye(3) -[10;20;30]];[0 0 0 1]];
    P2 = P1*Ht;
    He = H4xP2P(P1,P2);
    test.test = {'one camera -> exact & scale 1'};
    test.ok = max(max(abs(He\Ht-eye(4))))<1e-12;
    % two translated and scaled cameras
    C1 = {[10;10;10] [11;10;10]};
    R1 = {a2r([1;0;0],pi/10) a2r([1;0;0],-pi/10)};
    K1 = {[1000 0 500;0 1000 500;0 0 1] [1000 0 500;0 1000 500;0 0 1]};    
    C2 = {[0;0;0] [100;0;0]};
    R2 = R1;
    K2 = K1;
    P1 = cellfun(@(k,r,c) k*[r -r*c],K1,R1,C1,'UniformOutput',false);
    P2 = cellfun(@(k,r,c) k*[r -r*c],K2,R2,C2,'UniformOutput',false);
    He = H4xP2P(P1,P2,'fitCfitR');
    Ht = inv([100*[eye(3) [-10;-10;-10]];0 0 0 1]);
    test.test{end+1} = 'fitCfitR - two translated and scaled cameras';
    test.ok(end+1) = max(max(abs(He\Ht-eye(4))))<1e-12;
    He = H4xP2P(P1,P2,'setR&|b|');
    test.test{end+1} = 'setR&|b| - two translated and scaled cameras';
    test.ok(end+1) = max(max(abs(He\Ht-eye(4))))<1e-12;
    % two translated scaled camera centers and slightly rotated cameras around the baseline
    R1 = {a2r([1;0;0],pi/1000) a2r([1;0;0],-pi/1000)};
    R2 = {eye(3) eye(3)};
    P1 = cellfun(@(k,r,c) k*[r -r*c],K1,R1,C1,'UniformOutput',false);    
    P2 = cellfun(@(k,r,c) k*[r -r*c],K2,R2,C2,'UniformOutput',false);
    He = H4xP2P(P1,P2,'fitCfitR');
    test.test{end+1} = 'fitCfitR - two translated scaled camera centers and slightly rotated cameras around the baseline';
    test.ok(end+1) = max(max(abs(He\Ht-eye(4))))<1.572;
    H = test;
end

