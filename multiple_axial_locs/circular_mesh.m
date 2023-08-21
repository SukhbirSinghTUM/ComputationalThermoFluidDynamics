function [X , Y] = circular_mesh(M, r, h, theta)

dimR = size(M,1);
dim_theta = size(M,2);

R1 = r ; % inner radius 
R2 = r + h ;  % outer radius
nR = linspace(R2, R1, dimR) ;
nT = linspace(0,theta,dim_theta);
[R, T] = meshgrid(nR,nT) ;

X = R.*sind(T); 
Y = R.*cosd(T);


X = X';
Y= Y';

end