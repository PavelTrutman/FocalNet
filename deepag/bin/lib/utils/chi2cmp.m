%CHI2CMP  Pajdla: Compares a distribution of squared residuals with a chi2 distribution
%
%	function chi2cmp(res2,ni)
%
%	res	= column vector of squared residuals res2 = res.^2
%	ni	= degrres of freedom for chi^2 distribution
%
%       See also:  CHI2RND, QQPLOT, CUMHIST. 

%	Author: 	Tomas Pajdla, Tomas.Pajdla@esat.kuleuven.ac.be 
%					    pajdla@vision.felk.cvut.cz
%			11/02/94 ESAT-MI2, KU Leuven
%	Documentation:                 	 	  
%	Language: 	Matlab 4.2, (c) MathWorks
%       Last change  : $Id: chi2cmp.m,v 1.1 2005/04/28 16:54:38 pajdla Exp $
%       Status       : Ready ($Source: /home/cvs/Matlab/utils/chi2cmp.m,v $)  			 
%
function f = chi2cmp(r2,ni)

 len   = length(r2);
 hbins = max(10,len/20);

 if len > 1  

  stdr     = sqrt(sum(r2)/(len-1));
  sr2      = r2/stdr^2;

  meansr2  = mean(sr2);
  stdsr2   = std(sr2);
  maxsr2   = max(sr2);

  chi      = chi2rnd(ni*ones(size(sr2)));
  meanchi  = mean(chi);
  stdchi   = std(chi);

 else
  error('chi2cmp:Warning - data contain only 1 element');
 end	

 subplot(2,2,1);
 cumhist(sr2,hbins);
 axis([0 maxsr2 0 1]);
 axis('square');
 xlabel('squared standardized residuals');
 title('Histogram of an unknown data sample');

 subplot(2,2,2);
 cumhist(chi,hbins);
 axis([0 maxsr2 0 1]);
 axis('square'); 
 xlabel('squared standardized residuals');
 title('Histogram of the chi-square distribution');

 subplot(2,2,3)
 qqplot(sr2,chi,[10:10:100]);
 axis([0 Inf 0 Inf]);
 axis('square');
 xlabel('sres^2 Quantiles');
 ylabel('schi^2 Quantiles');
 title('Quantiles-Quantiles plot - [1:1:100] %');
 grid;

 subplot(2,2,4)
 dm = [ len;
        ni;
        stdr;
        meansr2;
        meanchi;
        stdsr2;
        stdchi];
 rlabel = ['data length        ';
           'degrees of freedom ';
           'std(res)           ';
           'mean(sres^2)       ';
           'mean(schi^2)       ';
           'std(sres^2)        ';
           'std(schi^2)        '];
 plottab(dm,rlabel,' ',' ','%.2f',12); 

return
