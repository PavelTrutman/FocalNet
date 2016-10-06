% Calculate feature vector - fundamental matrix - from correspondences
% size(u):=[2 number_of_points]
% (size(u1)==size(u2)):=true
%
% Oleh Rybkin, rybkiole@fel.cvut.cz
% INRIA, 2016

function [F,A] = F_features(u1, u2,method)
    if nargin<3
        method='Free';
    end    
    %some differents formats of input
    if size(u1,1)<3
       u1=[u1; ones(1,size(u1,2))];
       u2=[u2; ones(1,size(u1,2))];
    end
    %cannot calculate
    if size(u1,2)<7 || (size(u1,2)<8 && strcmp(method,'Free'))
        F=NaN;
        A=NaN;
        adprintf({},['they gave me ' num2str(size(u1,2)) ' points']);
        return;
    end
    % Fundamental matrix
    %[F,A]  = uu2F({u1,u2},{'None',method});
    [F,A]  = uu2F({u1,u2},{'[-1,1]',method});
    % choose if more then 1 answer
    for i=1:size(F,3) 
        F (:,:,i) = F(:,:,i)/norm(F(:,:,i));
        e(i) = max(abs(sum(u2.*(F(:,:,i)*u1))));
    end
    [~,ie] = min(e);
    F = F(:,:,ie);
    % normalization
    F = inv(A{2})'*F*inv(A{1});
    F=reshape(F,9,1);
    F = F/norm(F); [~,mi] = max(abs(F)); F = F*sign(F(mi));
end