function [estimatedI,xhat,r] = getSupport(A,K,y,en)
%[I,xhat,r] = getSupport(A,K,y)
%IN
%   A           Alphabet (M x dsize matrix)
%   K           Size of the support
%   y           Observation (vector of size M)
%OUT
%   I           Support
%   xhat        Estimate xhat= A*y + r
%   r           Residual

I= [];                               %Set of indices step after step
r= [];                               %set of residuals step after step

[M,dsize]= size(A);

xhat= zeros(dsize,1);
%Initialization
k=1;
[~,index]= sort(abs(A'*y),'descend');
I(:,k)= index(1:K);

%Convex problem 
cvx_begin quiet
    variable xhat_tmp(K,1)
    minimize( norm(y - A(:,I(:,k))*xhat_tmp )  )
    subject to
        A(:,I(:,k))*xhat_tmp >= eps
cvx_end

r(:,k)= y- A(:,I(:,k))*xhat_tmp; %(pinv(A(:,I(:,k)))*y);

while(1)
%    norm( r(:,k),2 )
   
   k=k+1;
   [~,index]= sort(abs(A'*r(:,k-1)),'descend');
   Ip= index(1:K);
   
   Iu= union(I(:,k-1),Ip);
   %Convex problem
   cvx_begin quiet
    variable xhat_tmp(length(Iu),1)
    minimize( norm(y - A(:,Iu)*xhat_tmp) )
    subject to
        A(:,Iu)*xhat_tmp >= eps
    cvx_end

   xhat(Iu,1)= xhat_tmp;%pinv(A(:,Iu))*y;
   xhat(setdiff(1:dsize,Iu))=0;
%    figure, plot(y); hold on; plot(A*xhat)

   [~,index]= sort(abs(xhat),'descend');
   I(:,k)= index(1:K);
   
   %Convex problem
   cvx_begin quiet
   variable xhat_tmp(K,1)
   minimize(  norm(y - A(:,I(:,k))*xhat_tmp ) )
   subject to
        A(:,I(:,k))*xhat_tmp >= eps
   cvx_end
   
   r(:,k)= y - A(:,I(:,k))*xhat_tmp;%(pinv(A(:,I(:,k)))*y);
%    figure, plot(y); hold on; plot(A(:,I(:,k))*xhat_tmp)
   


    if  ( norm(r(:,k)) <= norm(en) ) 
%         k= k - 1;
        estimatedI= I(:,k);
        xhat= zeros(dsize,1);
        xhat(estimatedI)= pinv(A(:,estimatedI))*y;
        break;
    elseif ( norm(r(:,k),2) >= norm(r(:,k-1),2) ) || k>2*K 
        k= k - 1;
        estimatedI= I(:,k);
        xhat= zeros(dsize,1);
        xhat(estimatedI)= pinv(A(:,estimatedI))*y;
        break;
    end
end
