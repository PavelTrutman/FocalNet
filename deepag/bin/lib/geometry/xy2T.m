% T = xy2T(x,y,tr) - 2D point registration transforms
%
% tr ... 'E' ~ euclidean 
%              x(:,1) -> y(:,1) & dir(x(:,1),x(:,2)) -> dir(y(:,1),y(:,2))
%    ... 'S' ~ similarity
%              x(:,1) -> y(:,1) & x(:,2) -> y(:,2) 
%    ... 'A' ~ affine
%              x(:,1) -> y(:,1) & x(:,2) -> y(:,2) & x(:,3) -> y(3,:)
%    ... 'H' ~ homography
%              x(:,1) -> y(:,1) & x(:,2) -> y(:,2) & x(:,3) -> y(:,3) & x(:,4) -> y(:,4)
%
% T  ... 3 x 3 matrix 
%
%        y = h2a(T*a2h(x))
%
% Run without parameters for demo.

% (c) T. Pajdla, pajdla@neovision.cz
% 13 Aug 2007
function T = xy2T(x,y,t)

if nargin > 0
    switch t
        case 'E' % euclidean
            x0 = x(:,1);
            y0 = y(:,1);
            vx = x(:,2)-x(:,1);
            vy = y(:,2)-y(:,1);
            R = xy2R(vx,vy);
            % y = R*(x-x0)+y0 = R*x + (y0-R*x0)
            T = [R                  y0-R*x0
                zeros(1,size(R,2))       1];
        case 'S' % similarity
            x0 = x(:,1);
            y0 = y(:,1);
            vx = x(:,2)-x(:,1); nx = norm(vx);
            vy = y(:,2)-y(:,1); ny = norm(vy);
            R = xy2R(vx,vy);
            % y = ny/nx*R*(x-x0)+y0 = ny/nx*R*x + (y0-ny/nx*R*x0)
            T = [ny/nx*R            y0-ny/nx*R*x0
                zeros(1,size(R,2))             1];
        case 'A' % affine
            T = a2h(y(:,1:3))*inv(a2h(x(:,1:3)));
        case 'H' % homography, see xy2H
            x = a2h(x(:,1:4));
            y = a2h(y(:,1:4));
            A = [[ zeros(size(y,2),3) -y([3 3 3],:)' y([2 2 2],:)'].*x([1 2 3 1 2 3 1 2 3],:)'
                [ y([3 3 3],:)' zeros(size(y,2),3) -y([1 1 1],:)'].*x([1 2 3 1 2 3 1 2 3],:)'
                [-y([2 2 2],:)' y([1 1 1],:)' zeros(size(y,2),3) ].*x([1 2 3 1 2 3 1 2 3],:)'];
            [V,D] = eig(A'*A); % faster but maybe less robust
            T = reshape(V(:,1),3,3)';
        otherwise
            error('xy2T: unknown transform type')
    end
else % demo
    x = rand(2,4);
    y = rand(2,4);
    %
    t = ['E';'S';'A';'H'];
    for i = 1:size(t,1)
        T = xy2T(x,y,t(i));
        z = h2a(T*a2h(x));
        %
        subfig(3,4,i);
        plot(x(1,:),x(2,:),'b.'); hold
        text(x(1,:)+0.01,x(2,:),num2str([1;2;3;4]),'color','b');
        plot(y(1,:),y(2,:),'r.','markersize',15);
        text(y(1,:)+0.02,y(2,:),num2str([1;2;3;4]),'color','r')
        plot(z(1,:),z(2,:),'g.');
        text(z(1,:)-0.02,z(2,:),num2str([1;2;3;4]),'color','g');
        plot(z(1,1:2),z(2,1:2),'g');
        axis equal
        title(t(i));
    end
end