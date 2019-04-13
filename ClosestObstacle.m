% Written by: Florence

function min_dis = ClosestObstacle(C,d,O)
%CLOSESTOBSTACLE returns closest obstacle 
%input:
%   C - ellipse parameter (nxn matrix)
%   d - ellipse parameter (n vector)
%   O - list of obstacles
%output:
%   min_dis - index of obstacle whose distance to the center of 
%   the ellipse is minimum 
Cinv = C^-1;

number_obs = size(O,3);


min_distance = zeros(1,number_obs);
for idx = 1:number_obs
    O(:,1,idx) = O(:,1,idx) - d(1);
    O(:,2,idx) = O(:,2,idx) - d(2);
    
    O(:,:,idx) = (Cinv* O(:,:,idx)')'; %transform into ball space
    min_distance(1,idx) = min(vecnorm(O(:,:,idx)')); %calculate the minimum distance
end

[d, min_dis] = min(min_distance);

end


