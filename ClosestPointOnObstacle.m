% Written by: Florence

function [ x ] = ClosestPointOnObstacle( C,d,O,i )
%CLOSESTPOINTONOBSTACLE returns closest point on the obstacle (create obstacle class)
%input:
%   C - ellipse parameter (nxn matrix)
%   d - ellipse parameter (n vector)
%   O - list of obstacles
%   i - index of obstacle closest to the ellipse
%output:
%   x - point on obstacle O that is closest to the ellipse

V=O(:,:,i)';
k=size(O,1); % dimension of vertices
n=size(O,2); % number of vertices

% transform each vertex to the ball space
V_bar=inv(C)*(V-d);

% define and solve problem in cvx
cvx_begin quiet
    variable x_bar(n) 
    variable w(k) nonnegative
    minimize(norm(x_bar))
    subject to
        V_bar*w==x_bar;
        sum(w)==1;
cvx_end

% transform point back to ellipse space
x=C'*x_bar+d;

end
