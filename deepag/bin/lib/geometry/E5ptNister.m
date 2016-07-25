% E = E5ptNister(X[,E]) - 5pt Minimal Relative Pose of Two Cameras by D. Nister 
%
% E        = essential matrix(ces) celarray
% X(1:2,:) = 3 x 5 points in image 1
% X(3:4,:) = 3 x 5 points in image 2

% T. Pajdla, pajdla@cmp.felk.cvut.cz
% 25 May 2015
function E = E5ptNister(X,E)

if nargin>0
    if nargin<2
        if all(size(X)==[6 5])
            % Exclude degeneracies due to repeated points
            s = [svd(X(1:2,:)) svd(X(4:5,:))]; s = s(end,:);
            if any(s<1e-8)
                E = [];
                return
            end
            % SOS to avoid failure for C1=C2 to get F with correct R and random t
            X = X+... 
                1e-8*[-0.0319    0.2822    0.2866   -0.1529    0.3569
                       0.2783   -0.4303   -0.0959   -0.3904   -0.0051
                       0.0000    0.0000    0.0000    0.0000    0.0000
                      -0.1352    0.2501    0.4205    0.3869   -0.2504
                      -0.3884   -0.3460    0.4984    0.4371    0.4825
                       0.0000    0.0000    0.0000    0.0000    0.0000]; 
            e  = solveE_nister(X(4:6,:),X(1:3,:)); % does not work for R ~ I
            eI = E4RIdentity(X(4:6,:),X(1:3,:)); % for R ~ I
            e(:,:,end+1) = eI;
            E = cell(size(e,3),1);
            for i=1:size(e,3)
                E{i} = e(:,:,i)/norm(e(:,:,i));
            end
        end
    end
else % unit tests
    E = [];
    % General scene
    X = [-2  2  2  2 -2
         -2 -2  2  2  2
          6  6  6  8  8];  
    i = 1; % test: General cameras  
    R1{i} = a2r([1;1;1], 1/10); R2{i} = a2r([1;1;1],-1/10);
    C1{i} = [-1;0;2]; C2{i} = [ 1;0;-5];
    i = 2; % test: P1 = [I|0] 
    R1{i} = a2r([1;1;1], 0); R2{i} = a2r([1;1;1],-1/10);
    C1{i} = [-1;0;2]; C2{i} = [ 1;0;-5];
    i = 3; % test: P2 = [I|0] 
    R1{i} = a2r([1;1;1], 1/10); R2{i} = a2r([1;1;1], 0);
    C1{i} = [-1;0;2]; C2{i} = [1;0;-5];    
    i = 4; % test: R2 = R1 
    R1{i} = a2r([1;1;1], 1/10); R2{i} = a2r([1;1;1], 1/10);
    C1{i} = [-1;0;2]; C2{i} = [1;0;-5];        
    i = 5; % test: C2 = C1
    R1{i} = a2r([1;1;1], 1/10); R2{i} = a2r([1;1;1], -1/10);
    C1{i} = [-1;0;2]; C2{i} = [-1;0;2];            
    i = 6; % test: C1 = C2 & P1 = [ I | 0] 
    R1{i} = a2r([1;1;1], 0); R2{i} = a2r([1;1;1],-1/10);
    C1{i} = [0;0;0]; C2{i} = [0;0;0];    
    % Run the test
    for i=1:numel(R1)
        P1 = [R1{i} -R1{i}*C1{i}]; P2 = [R2{i} -R2{i}*C2{i}];
        x = [a2h(h2a(P1*a2h(X)));a2h(h2a(P2*a2h(X)))];
        Ft = PP2F(P1,P2); Ft = Ft/norm(Ft); % normalized F contructed from the projection matrices
        F = E5ptNister(x); % estimated F
        if ~isempty(F)
            d = min(cellfun(@(x) max(min(abs([x(:)-Ft(:) x(:)+Ft(:)]),[],2)),F)); % matrix distance by the largest element absolute difference considering unknown sign of F
            E(1,i) = d<1e-6; % is OK?
            p2 = cellfun(@(x) x{2},cellfun(@(x) E2PP(x),F,'UniformOutput',false),'UniformOutput',false); p2 = cat(2,p2{:}); % P2 matrices
            r2 = cellfun(@(x) P2R(x),p2,'UniformOutput',false); %rotations
        else
            E(1:2,i) = false;
        end
    end
end

function e = E4RIdentity(y,x)
%A = [cr(x2(:,1),x1(:,1)) cr(x2(:,2),x1(:,2)) cr(x2(:,3),x1(:,3)) cr(x2(:,4),x1(:,4)) cr(x2(:,5),x1(:,5))]';
 A = [x(2,1)*y(3,1)-x(3,1)*y(2,1) -x(1,1)*y(3,1)+x(3,1)*y(1,1) x(1,1)*y(2,1)-x(2,1)*y(1,1)
      x(2,2)*y(3,2)-x(3,2)*y(2,2) -x(1,2)*y(3,2)+x(3,2)*y(1,2) x(1,2)*y(2,2)-x(2,2)*y(1,2)
      x(2,3)*y(3,3)-x(3,3)*y(2,3) -x(1,3)*y(3,3)+x(3,3)*y(1,3) x(1,3)*y(2,3)-x(2,3)*y(1,3)
      x(2,4)*y(3,4)-x(3,4)*y(2,4) -x(1,4)*y(3,4)+x(3,4)*y(1,4) x(1,4)*y(2,4)-x(2,4)*y(1,4)
      x(2,5)*y(3,5)-x(3,5)*y(2,5) -x(1,5)*y(3,5)+x(3,5)*y(1,5) x(1,5)*y(2,5)-x(2,5)*y(1,5)];
 %[~,~,e] = svd(A);
 [e,~] = eig(A'*A);
 e = [0 -e(3,1) e(2,1);e(3,1) 0 -e(1,1);-e(2,1) e(1,1) 0];
return
function z = cr(x,y)
  z = [x(2)*y(3)-x(3)*y(2);-x(1)*y(3)+x(3)*y(1);x(1)*y(2)-x(2)*y(1)];
return
function R = P2R(P)
  [K,R] = P2KRC(P);
return
