% Calculate feature vector from correspondences
% 
% used in deepag, change this function to change format of deepag.m output
% in:
%     x - matrix of 2d coordinates in the first image of the form size(x):=[number_of_points 2]
%     p - matrix of coordinates in the second image corresponding to x
% out:
%     Mvec - feature vector
%
% Oleh Rybkin, rybkiole@fel.cvut.cz
% INRIA, 2016

function [Mvec] = get_features(x, xp)

    %placeholder
    if 1
        Mvec=1;
        return;
    end

    %use Fundamental matrix formulation
    if 0
        [Mvec,~]=F_features(x',xp');
        return;
    end
    
    %use formulation via coefficients to monomials
    if 0
        Mvec=monomial_features(x,xp)
    end;
end