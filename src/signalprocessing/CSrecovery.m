clear all
close all
clc

%Choose problem parameters
M=26;
dsize=50;
K=5;
Ampl_max=20;

%Create data
A=randn(M,dsize)*sqrt(1/M)+2;         %Dictionary
A=A./(ones(M,1)*sqrt(sum(A.^2)));   %Normalization
x=zeros(dsize,1);                   %Sparse vector

trueI= round((dsize-1)*rand(K,1))+1;%nz indices
x(trueI)= Ampl_max*(rand(K,1)-0); %nz values
y= A*x;                             %Observation

%% Now we want to find the sparse vector x from y, K and A
[I,xhat,r] = getSupport(A,K,y);
figure, plot(y); hold on; plot(A*xhat);
legend('y','yhat');


figure, stem(x,'LineWidth',2); hold on; stem(xhat);
legend('x','xhat');