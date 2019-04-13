% Written by: Florence

function [ a,b ] = tangent( C,d,x )
%tangentPlane Calculate a hyperplane tangent to the ellipse defined by
%               C, d at x.
% input:
%   C - symmetric, semi-positive definite matrix
%   d - n-vector, displacement of ellipse
%   x - n-vector, point on obstacle closest to ellipse
% output:
%   a - n-vector, slope of hyperplane
%   b - n-vector, root of hyperplane

a=2*inv(C)*inv(C')*(x-d);
b=a.'*x;
end
