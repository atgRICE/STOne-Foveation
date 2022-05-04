%  Function createDifferenceOperators by Tom Goldstein
%           modified to createDifferenceOperators_fov by Anthony Giljum
%
%  This function creates first-order difference operators
%  in the x- and y-direction.  The operators act of images that  
%  have been vectorized using the nested dissection ordering.  The 
%  returned values are sparse matrices that perform the differencing.

function [ Dx, Dy ] = createDifferenceOperators_fov(order,order_prev,topLeft,botRight)

rows = order.n;
cols = order.n;

%  Indexes to the location of each pixel in the vector
index_full = order.matrix*(order_prev.n/order.n)^2;
index_prev = kron(order_prev.matrix,ones(order.n/order_prev.n));
        % We now have that ceil(index_full)==index_prev
index_tmp = index_prev;
for i = 1:size(topLeft,1)
    index_tmp(topLeft(i,1):botRight(i,1),topLeft(i,2):botRight(i,2)) = ...
                index_full(topLeft(i,1):botRight(i,1),topLeft(i,2):botRight(i,2));
end
index = index_tmp*(order.n/order_prev.n)^2;

% Label each block with an integer, independent of its resolution, in
% accordance with the nested embedding.
tmp = index(:);
tmpsort = sort(unique(tmp))';
lookuptable(tmpsort) = 1:length(tmpsort);
tmp = lookuptable(tmp);
index = reshape(tmp,size(index));

% Build the derivative operator in the x (rows) direction
xInd = index; %index of each pixel
diffInd = circshift(index,[-1 0]); % index of adjacent pixels
vals = ones(rows,cols);
vals(end,:)=0;  % don't take differences along the edges
valsB = vals;
for i = 1:size(botRight,1)
    vals(topLeft(i,1):botRight(i,1),topLeft(i,2):botRight(i,2)) = order.n/order_prev.n;
    valsB(botRight(i,1),topLeft(i,2):botRight(i,2)) = order_prev.n/order.n;
    valsB(topLeft(i,1):botRight(i,1),topLeft(i,2):botRight(i,2)) = order.n/order_prev.n;
end
xInd = xInd(:);
diffInd=diffInd(:);
vals = vals(:);
valsB = valsB(:);
Dx = sparse([xInd;xInd],[xInd;diffInd], [-vals;valsB]); % Build sparse matrix
Dx = Dx / (order.n/order_prev.n);

% Build derivative operator in the y (cols) direction
yInd = index;
diffInd = circshift(index,[0 -1]);
vals = ones(rows,cols);
vals(:,end)=0;  % don't take differences along the edges
valsR = vals;
for i = 1:size(botRight,1)
    vals(topLeft(i,1):botRight(i,1),topLeft(i,2):botRight(i,2)) = order.n/order_prev.n;
    valsR(topLeft(i,1):botRight(i,1),botRight(i,2)) = order_prev.n/order.n; % RHS Border
    valsR(topLeft(i,1):botRight(i,1),topLeft(i,2):botRight(i,2)) = order.n/order_prev.n;
end
yInd = yInd(:);
diffInd=diffInd(:);
vals = vals(:);
valsR = valsR(:);
Dy = sparse([yInd;yInd],[yInd;diffInd], [-vals;valsR]); % Build sparse matrix
Dy = Dy / (order.n/order_prev.n);

return;
