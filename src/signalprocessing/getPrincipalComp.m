function [vector_, index ] = getPrincipalComp(vector, nbprincipal)
%function [ lowdim_index ] = getPrincipal(nbprincipal,zhat)
%   INPUT
%         nbprincipal:   integer between 1 and length(vector)
%         vector:        line or column
%   OUTPUT
%         vector_:       initial vector with only nbprincipal nz components
%         index:         index of the sorted initial vector

[~,index]=sort(vector.^2,'descend');
vector_=zeros(length(vector),1);

vector_(index(1:nbprincipal))= vector(index(1:nbprincipal));

end

