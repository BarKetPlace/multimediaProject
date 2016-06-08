clear all
close all
clc

%Choose problem parameters
M=26;
dsize=50;
K=25;
Ampl_max=20;

%Create data
A=randn(M,dsize)*sqrt(1/M)+2;         %Dictionary
A=A./(ones(M,1)*sqrt(sum(A.^2)));   %Normalization
x=zeros(dsize,1);                   %Sparse vector

trueI= round((dsize-1)*rand(K,1))+1;%nz indices
x(trueI)= Ampl_max*(rand(K,1)-0); %nz values
y= A*x;                             %Observation

%% Now we want to find the sparse vector x from y, K and A
I= [];                               %Set of indices step after step
r= [];                               %set of residuals step after step

xhat= zeros(dsize,1);
%Initialization
k=1;
[~,index]= sort(A'*y,'descend');
I(:,k)= index(1:K);

%Convex problem 
cvx_begin quiet
    variable xhat_tmp(K,1)
    minimize( norm(y - A(:,I(:,k))*xhat_tmp )  )
    subject to
        xhat_tmp >= eps
cvx_end

r(:,k)= y- A(:,I(:,k))*xhat_tmp %(pinv(A(:,I(:,k)))*y);
% figure, plot(y); hold on; plot(A(:,I(:,k))*xhat_tmp)
while(1)
   norm( r(:,k),2 )
   
   k=k+1;
   [~,index]= sort(A'*r(:,k-1),'descend');
   Ip= index(1:K);
   
   Iu= union(I(:,k-1),Ip);
   %Convex problem
   cvx_begin quiet
    variable xhat_tmp(length(Iu),1)
    minimize( norm(y - A(:,Iu)*xhat_tmp) )
    subject to
        xhat_tmp >= 0
    cvx_end

   xhat(Iu,1)= xhat_tmp;%pinv(A(:,Iu))*y;
   xhat(setdiff(1:dsize,Iu))=0;
%    figure, plot(y); hold on; plot(A*xhat)

   [~,index]= sort(xhat,'descend');
   I(:,k)= index(1:K);
   
   %Convex problem
   cvx_begin quiet
   variable xhat_tmp(K,1)
   minimize( norm(y - A(:,I(:,k))*xhat_tmp )  )
   subject to
        xhat_tmp >= eps
   cvx_end
   
   r(:,k)= y - A(:,I(:,k))*xhat_tmp;%(pinv(A(:,I(:,k)))*y);
%    figure, plot(y); hold on; plot(A(:,I(:,k))*xhat_tmp)
   


    if  ( norm(r(:,k),2) > norm(r(:,k-1),2) ) || k>2*K
        k= k - 1;
        estimatedI= I(:,k);
        xhat= zeros(dsize,1);
        xhat(estimatedI)= pinv(A(:,estimatedI))*y;
        break;
    end
end

figure, plot(y); hold on; plot(A*xhat);
legend('y','yhat');


figure, stem(x,'LineWidth',2); hold on; stem(xhat);
legend('x','xhat');