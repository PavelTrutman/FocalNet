% a = PPX2ae(P1,P2,X) - Apical angles of the epipolar triangles P1, X, P2 
%
% P1, P2 = 3 x 4 camera projection matrix
% X = 3 x n 3D points
% a = apical angles [rad] in the epipolar triangle C1-X-C2
%     a(1,:) = angle(C2,C1,X) ... apical angle at camera 1
%     a(2,:) = angle(X,C2,C1) ... apical angle at camera 2
%     a(3,:) = angle(C1,X,C2) ... apical angle at X

% T. Pajdla, pajdla@cmp.felk.cvut.cz, 2015-09-05
function a = PPX2ae(P1,P2,X)
if nargin>0
    if ~isempty(X)
        %[~,~,C1] = P2KRC(P1); % camera center 1 ... slow
        %[~,~,C2] = P2KRC(P2); % camera center 2
        C1 = -P1(:,1:3)\P1(:,4);
        C2 = -P2(:,1:3)\P2(:,4);
        v1 = C1-C2;  v1 = v1/vnorm(v1); % C2->C1        
        C1 = repmat(C1,1,size(X,2));
        C2 = repmat(C2,1,size(X,2));
        v2 = X-C1; v2 = v2./repmat(vnorm(v2),[size(v2,1) 1]); % C1->X
        v3 = C2-X; v3 = v3./repmat(vnorm(v3),[size(v3,1) 1]); % X->C2
        a = [v1'*v2;v1'*v3];
        a = pi-acos(a);
        a(3,:) = pi-sum(a);
    else
        a = [];
    end
else % unit tests
    test.fname = mfilename;
    test.test = {'Known values'};
    P1 = [eye(3) -[0;0;0]];
    P2 = [eye(3) -[2;0;0]];
    X  = [0 1 2 1   1   1
          1 1 1 1/2 1/3 4
          0 0 0 0   0   0];
    a = PPX2ae(P1,P2,X);    
    test.ok = all(abs(a(3,:)-[ 1.1071 1.5708 1.1071 2.2143 2.4981 0.4900])<1e-4);
    a = test;
end