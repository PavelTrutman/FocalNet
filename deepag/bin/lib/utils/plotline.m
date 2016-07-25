% h = plotline(l[,rect]) - Line plotting clipped by a rectangle
%
% l         ... 3 x n line vectors
% rect      ... clipping rectangle [minx maxx miny maxy]
%               obtained from axis if missing, uses the first two coordinates
% h         ... line hadles
%
% See also LINE, AXIS

% (c) T.Pajdla, www.neovision.cz, Sep 17, 2005
function [h vx] = plotline(l,r)

% prevent division by zero warnings but save the curent warning state
warn = warning('off','MATLAB:divideByZero');

if nargin<2
    r = axis;
    r = r(1:4);
end
% normalized line vectors
l = l./([1;1;1]*sqrt(sum(l.^2)));
% normalized rectangle vertices
rv = a2h(r([1 3;2 3;2 4;1 4])');
rv = rv./([1;1;1]*sqrt(sum(rv.^2)));
% normalized rectagle lines
rl = [cross(rv(:,1),rv(:,2)), cross(rv(:,2),rv(:,3)), cross(rv(:,3),rv(:,4)), cross(rv(:,4),rv(:,1))];
rl = rl./([1;1;1]*sqrt(sum(rl.^2)));
% the center of the rectangle
c  = mean(reshape(r,2,2))';
dm = sqrt(sum((c-h2a(rv(:,1))).^2)); % maximal distance
k  = 1;
vx = [];
for i=1:size(l,2)
    % intersection vertices
    v = [cross(l(:,i),rl(:,1)), cross(l(:,i),rl(:,2)), cross(l(:,i),rl(:,3)), cross(l(:,i),rl(:,4))];
    v = h2a(v);
    % pairwise distances
    vd = sqrt(sum((kron(v,[1 1 1 1])-kron([1 1 1 1],v)).^2));
    vd = reshape(vd,4,4);
    % symmetric => keep only the values above the main diagonal
    vd(logical(tril(ones(4)))) = inf;
    % the close points
    ivd = vd<100*eps;
    % remove points that are too close to some other points
    rmi = [];
    for j=2:size(ivd,2)
        if any(ivd(:,j)) 
            rmi = [rmi ; j];
            ivd(j,:) = false;
        end
    end
    kpi = true(1,size(ivd,2)); % keep those that should not be removed
    kpi(rmi) = false;
    v = v(:,kpi);
    % the distances from the center
    d = sqrt(sum([v-c*ones(1,size(v,2))].^2));
    % keep only the ones close to the center
    dix = d<(dm+100*eps);
    d = d(dix);
    v = v(:,dix);
    if length(d)>=2
        % orderd by the distance from the center
        [sd,si] = sort(d); 
        % choose the two closest ones
        vx(:,:,k) = v(:,si(1:2)); 
        k         = k + 1;
    end
end
if ~isempty(vx)
    h = line(squeeze(vx(1,:,:)),squeeze(vx(2,:,:)));
else
    h = [];
end
% restore the waring state
warning(warn.state,warn.identifier);
