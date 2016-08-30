% [F,A] = uu2F(u[,mth]) - Fundamental matrix computation
%
% u   = {u1 u2}, image matches u2'*F*u1 = 0
% mth = method ({'HZ','Free'} implicit)
%       mth{1} = normalization
%                'None' = no normalization
%                'HZ' = full affine as in HZ-2003
%                '[-1,1]' = ranges of points are scaled to interval [-1,1]x[-1,1]
%       mth{2} = constraints on F imposed
%                'Free' = no constraints, i.e. |F|=0 not required
% F   = Fundamental matrix u2'*F*u1 = 0
% A   = normalization transforms: F = A{2}'*Fn*A{1};  u{2}'*A{2}'*Fn*A{1}*u{1}->0 => Fn
%
% T. Pajdla, pajdla@cvut.cz, 2016-08-28
function [F,A] = uu2F(u,mth)
if nargin>0
    if nargin<2
        mth = {'HZ','Free'};
    end
    if ~iscell(mth)
        mth = {mth 'Free'};
    end
    u = cellfunuo(@(x) x(1:2,:),u); % use the first two coordinates
    switch mth{1}
        case 'None'
            A{1} = eye(3); A{2} = A{1};
        case 'HZ'
            [A{1},A{2}] = xy2nxy(u{1},u{2});
        case '[-1,1]'
            A{1} = x2nx(u{1},'[-1,1]');
            A{2} = x2nx(u{2},'[-1,1]');
    end
    % normalize
    switch mth{2}
        case 'Free'
            x = {A{1}*a2h(u{1}) A{2}*a2h(u{2})};
            B = zeros(size(x{1},2),9);
            % B * f = 0
            for i=1:size(x{1},2)
                B(i,:) = m2v(x{1}(:,i)*x{2}(:,i)');
            end
            % f is the null space of B
            [~,~,f] = svd(B,0);
            f = f(:,end);
            % reformat back
            F = reshape(f,3,3)';
            % denormalize: u2'*Fa*u1 = u2'*A2'*F*A1*u1 => Fa = A2'*F*A
            F = A{2}'*F*A{1};
        otherwise
            error([mth{2} ' not implemented']);
    end
else % unit tests
    % test 1 
    X = [0 1 1 0 0 1 2 0
         0 0 1 1 0 0 1 1
         0 0 0 0 1 2 1 2];
    P1 = [1 0 0 1
          0 1 0 0
          0 0 1 1];
    P2 = [1 0 0 0
          0 1 0 1
          0 0 1 1];
    u1 = X2u(X,P1);
    u2 = X2u(X,P2);
    Fo = uu2F({u1,u2},{'[-1,1]','Free'}); 
    Fo = Fo/norm(Fo);
    F(1) = max(abs(sum(u2.*(Fo*u1))))<1e-8;
    % test 2
    E = E5ptNister([u1(:,1:5);u2(:,1:5)]);     
    F(2) = norm(E{end}-Fo/norm(Fo,2),2)<1e-8;
    % test 3
    P1 = rand(3,4);
    P2 = rand(3,4);
    u1 = X2u(X,P1);
    u2 = X2u(X,P2);
    Fo = uu2F({u1,u2},{'[-1,1]','Free'});
    Fo = Fo/norm(Fo);    
    F(3) = max(abs(sum(u2.*(Fo*u1))));
    % test 4
    P1 = [1 0 0 0
          0 1 0 0
          0 0 1 0];
    P1(:,4) = rand(3,1);  
    P2 = [1 0 0 0
          0 1 0 0
          0 0 1 0];  
    P2(:,4) = rand(3,1);
    X = rand(3,8);
    u1 = X2u(X,P1);
    u2 = X2u(X,P2);    
    Fo = uu2F({u1,u2},{'HZ','Free'});
    Fo = Fo/norm(Fo);
    E = E5ptNister([u1(:,1:5);u2(:,1:5)]);    
    F(4) = norm(E{end}/E{end}(1,3)-Fo/Fo(1,3),2)<1e8;
end

