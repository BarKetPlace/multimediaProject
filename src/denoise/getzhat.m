function [zhat] = getzhat(D,ey)
%function [zhat] = getzhat(D,ey)
%IN ::  D      -> Dictionnary
%       ey      -> Observation
%OUT::  zhat   -> estimate coefficients.
%The actual estimate of y is D*zhat

[M,dsize]=size(D); % Number of mfcc
SparsityTarget = .5; %
%Sparsity constrain
[B] = lasso(D,ey);
%Percentage of 0 coef
perc=sum(B==0)/dsize;
[~,idx]= min(abs(perc-SparsityTarget));
zhat= B(:,idx);

Dtrunc=D(:,zhat~=0);

%Positivity constrain
%We want to minimize ||ey-Dtrunc*zhat_lsq||^2 s.t. -Dtrunc*zhat_lsq <= 0
C=Dtrunc;
d=ey;
A=-Dtrunc;
b=zeros(1,M);
%No bound
ub=inf(size(C,2),1);
lb=-1*ub;

options = optimoptions('lsqlin','Algorithm','active-set');
[zhat_lsq]= lsqlin(C,d,A,b,[],[],lb,ub,[],options);

% eyhat=Dtrunc*zhat_lsq;
% tmp = find(eyhat<=0);
% if any(tmp)
%     warning('Non positive values');
% end
zhat(find(zhat~=0))=zhat_lsq;
end

