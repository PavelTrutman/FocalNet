% [F,A] = uu2F(u[,mth]) - Fundamental matrix computation
%
% u   = {u1 u2}, image matches u2'*F*u1 = 0
% mth = method ({'HZ','Free'} implicit)
%       mth{1} = normalization
%                'None' = no normalization
%                'HZ' = full affine as in HZ-2003
%                '[-1,1]' = ranges of points are scaled to interval [-1,1]x[-1,1]
%       mth{2} = constraints on F imposed
%                'Free'  = no constraints, i.e. |F|=0 not required (works for 8 and more matches)
%                '|F|=0' = rank two F imposed (works for 7 and more matches) 
% F   = Fundamental matrix u2'*F*u1 = 0
% A   = normalization transforms: F = A{2}'*Fn*A{1};  u{2}'*A{2}'*Fn*A{1}*u{1}->0 => Fn
%
% T. Pajdla, pajdla@cvut.cz, 2016-09-11
function [F,A] = uu2F(u,mth)
if nargin>0
    if nargin<2
        mth = {'HZ','Free'};
    end
    if ~iscell(mth)
        mth = {mth 'Free'};
    end
    if numel(mth)<2
        error('uu2F: cellarray mth must have two elements');
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
        otherwise
            error('uu2F: unknown normalization method %s',mth{1});
    end
    % normalize
    x = {A{1}*a2h(u{1}) A{2}*a2h(u{2})};
    B = zeros(size(x{1},2),9);
    % B * f = 0
    for i=1:size(x{1},2)
        B(i,:) = m2v(x{1}(:,i)*x{2}(:,i)');
    end
    % f is the null space of B
    [~,~,f] = svd(B,0);
    switch mth{2}
        case 'Free'
            if size(B,1)<8
                error('uu2F: mth{2} = ''Free'' needs at least 8 matches');
            end
            f = f(:,end);
            % reformat back
            F = reshape(f,3,3)';
        case '|F|=0'
            % det(F1+t*F2) = a'*[t^3;t^2;t;1]
            % det(t*F1+F2) = a'*[1;t;t^2;t^3]
            % Choose the one with the larger absolute value at the highest power
            % to avoid division by zero and avoid t going to infty
            N{1} = reshape(f(:,end-1),3,3)';
            N{2} = reshape(f(:,end),3,3)';
            a = [det(N{2}) 
                (det([N{1}(:,1) N{2}(:,[2 3])])-det([N{1}(:,2) N{2}(:,[1 3])])+det([N{1}(:,3) N{2}(:,[1 2])]))
                (det([N{2}(:,1) N{1}(:,[2 3])])-det([N{2}(:,2) N{1}(:,[1 3])])+det([N{2}(:,3) N{1}(:,[1 2])]))
                 det(N{1})];
            if abs(a(end))>abs(a(1)) % choose the larger value
                a = a(end:-1:1);
                ix = [1 2]; % det(t*F1+F2) = a'*[1;t;t^2;t^3]
            else
                ix = [2 1]; % det(F1+t*F2) = a'*[t^3;t^2;t;1]
            end
            a = a/a(1); % monic polynomial
            C = [0 0 -a(4);1 0 -a(3);0 1 -a(2)]; % companion matrix
            t = eig(C);
            t = t(abs(imag(t))<eps); % select real solutions
            F = zeros(3,3,numel(t));
            for i=1:numel(t)
                F(:,:,i) = t(i)*N{ix(1)}+N{ix(2)};
            end
        otherwise
            error([mth{2} ' not implemented']);
    end
    % denormalize: u2'*Fa*u1 = u2'*A2'*F*A1*u1 => Fa = A2'*F*A
    for i=1:size(F,3)
        F(:,:,i) = A{2}'*F(:,:,i)*A{1};
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
    Fo = uu2F({u1(:,1:8),u2(:,1:8)},{'[-1,1]','|F|=0'}); 
    for i=1:size(Fo,3)
        Fo(:,:,i) = Fo(:,:,i)/norm(Fo(:,:,i));
        e(i) = max(abs(sum(u2.*(Fo(:,:,i)*u1))));
    end
    F(3) = any(e<1e-10);    
    % test 4
    P1 = rand(3,4);
    P2 = rand(3,4);
    u1 = X2u(X,P1);
    u2 = X2u(X,P2);
    Fo = uu2F({u1,u2},{'[-1,1]','Free'});
    Fo = Fo/norm(Fo);    
    F(4) = max(abs(sum(u2.*(Fo*u1))))<1e-10;
    % test 5
    Fo = uu2F({u1,u2},{'[-1,1]','|F|=0'});
    for i=1:size(Fo,3)
        Fo(:,:,i) = Fo(:,:,i)/norm(Fo(:,:,i));
        e(i) = max(abs(sum(u2.*(Fo(:,:,i)*u1))));
    end    
    F(5) = min(e)<1e-10;    
    % test 6 & 7
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
    Fd = uu2F({u1,u2},{'[-1,1]','|F|=0'});
    for i=1:size(Fd,3)
        Fd(:,:,i) = Fd(:,:,i)/norm(Fd(:,:,i));
        e(i) = max(abs(sum(u2.*(Fd(:,:,i)*u1))));
    end    
    [~,ie] = min(e);
    Fd = Fd(:,:,ie);
    E = E5ptNister([u1(:,1:5);u2(:,1:5)]);    
    F(6) = norm(E{end}/E{end}(1,3)-Fo/Fo(1,3),2)<1e8;    
    F(7) = (norm(Fd/Fd(1,3)-Fo/Fo(1,3),2)<1e8) && (min(svd(Fd))<1e-14);    
end
