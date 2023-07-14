function [X,Y] = setUpMesh(M, l, formfunction)

dx = l/(size(M,2)-1);
X = repmat(0:dx:l,size(M,1),1);

Y = M;
Y(1,:) = formfunction(0:dx/l:1);
for i=1:size(M,2)
    Y(:,i) = (Y(1,i):-Y(1,i)/((size(M,1)-1)):0)';
end