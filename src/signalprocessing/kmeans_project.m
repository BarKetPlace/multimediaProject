function [idxs, C] = kmeans_project(DATA, num_clust )
%KMEANS_PROJECT Perform k-means clustering on DATA
%   INPUT  :    DATA      : nxm dimensional data. n = variables, m = number of
%                           points
%               num_clust : number of cluster
%   OUTPUT :    idx       : indices of cluster for each point
%               C         : centroids matrix
%   The function works fine if we compare to the built-in matlab function.
%   However it is slower, especially for high dimensionnal data.

    % Initialization of parameters
row     = size(DATA, 1) ;
column  = size(DATA, 2) ;
Nb_iter = 100 ;

    % 1) Initialization of centroids.
    % METHOD : choose randomly within the data (Forgy method)
C = DATA(:, randperm(column, num_clust)) ;          % get num_clust unique indices from all the possible indices

DATAmat = repmat(DATA, [1 1 num_clust]) ;           % initialize DATAmat matrix which will be used to compute distance
DATAmat = permute(DATAmat, [1 3 2]) ;               % permute dimension so that matrices fit

for j = 1 : Nb_iter
    % 2) For each point in the DATA, find indice of nearest centroid using euclidian
    % distance
Cmat = repmat(C, [1 1 column]) ;                    % add dimension to matrix to compute everything without for loop

dist = squeeze(sum((Cmat - DATAmat).^2)) ;          % compute euclidian distance between each point and each centroids

idxs = repmat(min(dist), [num_clust 1]) ;           % get the indice of the nearest centroid for each point
idxs = (1:num_clust) * (dist == idxs) ;             %

    % 3) Update the centroid of each cluster by computing the mean of all
    % the points contained in each cluster
for i = 1 : num_clust
   x = DATA(repmat(idxs == i, [row 1])) ;           % get all the points in cluster i
   x = reshape(x, [row length(x)/row]) ;            %
   
   C(:, i) = mean(x, 2) ;                           % computes the mean and stores it in centroid vector
end
end


end

