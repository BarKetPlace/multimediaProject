function [vector_, index,nbprincipal ] = getPrincipalComp(vector, energyperc)
%[vector_, index,nbprincipal ] = getPrincipalComp(vector, energyperc)
%   INPUT
%       vector:        line or column  
%       energyperc:    percentage of energy 
%         
%   OUTPUT
%         vector_:       initial vector with only nbprincipal nz components
%         index:         index of the sorted initial vector

[~,index]=sort(vector.^2,'descend');
vector_=zeros(length(vector),1);
totenergy= sum(vector.^2);
cumenergy= cumsum(vector(index).^2)./totenergy;
%The diffenergy is typicaly negative in the begining and then turns
%negative. The change in sign means that the energy becomes greater than
%the threshold
diffenergy = cumenergy-energyperc;

nbprincipal=1;
% If there is a change of sign in the difference vector
if (any(diffenergy<=0) && any(diffenergy)>=0)
    nbprincipal= find(diffenergy(1:end-1,1).*diffenergy(2:end,1)<0);
end 
    
vector_(index(1:nbprincipal))= vector(index(1:nbprincipal));

end

