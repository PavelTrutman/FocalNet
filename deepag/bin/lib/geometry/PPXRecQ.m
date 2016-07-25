% q = PPXRecQ(P1,P2,X,type,par) - Quality of the Reconstruction of X from cameras P1, P2, the bigger the better 
%
% P1, P2 = 3 x 4 camera projection matrix
% X      = 3 x n 3D points
% type   = type of the quality measure
%          'ApXMeanN' = number of points N times the mean of the apical angle at points X (implicit)
%          par = best angle (pi/4 if missing)
%          'ApXMediN' = number of points N times the median of the apical angle at points X
%          par = best angle (pi/4 if missing)
%          'ApX>parN' = the number of apical angles at X larger tan par(2)
%          par(1) = best angle (pi/4 if missing)   
%          par(2) = only apical angles >= par(2) are counted (pi/180 if missing)

% T. Pajdla, pajdla@cmp.felk.cvut.cz, 2015-09-05
function q = PPXRecQ(P1,P2,X,type,par)
if nargin<4
    type = 'ApXMeanN';
end
switch type(1:3)
    case 'ApX'
        if nargin<5
            par(1) = pi/4; % best angle
        end
        q = PPX2ae(P1,P2,X); % the apical angles
        if ~isempty(q)
            q = q(3,:); % at points X
            q(q>par(1)) = par(1)-q(q>par(1)); % par is the best angle
            switch type(4:7)
                case 'Mean'
                    q = length(q)*mean(q);
                case 'Medi'
                    q = length(q)*median(q);
                case '>par'
                    if length(par)<2
                        par(2) = pi/180;
                    end
                    q = sum(q>par(2));
                otherwise
                    error('unknown type of qulity measure');
            end
        else
            q = -Inf;
        end
    otherwise
        error('unknown type of qulity measure');
end


