function [zhat] = getzhat(D,ey)
%function [zhat] = getzhat(D,ey)
%IN ::  D      -> Dictionnary
%       ey      -> Observation
%OUT::  zhat   -> estimate coefficients.
%The actual estimate of y is D*zhat
numiter=1;
[M dsize]=size(D);
W = ones(dsize,1); % initial weights
delta= 10^-5;
lambda=.001;
for k=1:numiter
    
cvx_begin quiet
    variable zhat(dsize)
    
    minimize( norm( D * zhat - ey, 2 ) + lambda*norm( zhat, 1 ) )
    subject to
        D * zhat >= eps
cvx_end


%   W = 1./(delta + abs(zhat));
%       figure(1), clf;
%     subplot(121);
%     plot(ey,'LineWidth',2); hold on; plot(D*zhat);
%     subplot(122);
%     stem(zhat);
%     title(['nzz:: ' num2str(length(find( abs(zhat) > delta )))]);
end

end








% [M,dsize]=size(D); % Number of mfcc
% SparsityTarget = .5; %
% %Sparsity constrain
% [B] = lasso(D,ey);
% %Percentage of 0 coef in z
% perc=sum(B==0)/dsize;
% [~,idx]= min(abs(perc-SparsityTarget));
% zhat= B(:,idx);
% 
% Dtrunc=D(:,zhat~=0);
% 
% %Positivity constrain
% %We want to minimize ||ey-Dtrunc*zhat_lsq||^2 s.t. -Dtrunc*zhat_lsq <= -eps
% C=Dtrunc;
% d=ey;
% A=-Dtrunc;
% b=zeros(1,M)-eps;
% %No bound
% ub=inf(size(C,2),1);
% lb=-1*ub;
% %We choose the right algorithm that works for underdetermined system (less
% %equations than variable)
% options = optimoptions('lsqlin','Algorithm','active-set');
% 
% [zhat_lsq]= lsqlin(C,d,A,b,[],[],lb,ub,[],options);
% 
% % eyhat=Dtrunc*zhat_lsq;
% % tmp = find(eyhat<=0);
% % if any(tmp)
% %     warning('Non positive values');
% % end
% zhat(find(zhat~=0))=zhat_lsq;



