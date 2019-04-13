% Written by: Florence
% Run this to generate the test data

% Specify boundaries and number of obstacles
n=80;       % number of obstacles
limit=7.5;    % abs max of boundary

% create a template triangle
L = linspace(0,2*pi,4);
x = cos(L)';
y = sin(L)';

% apply random transformations to template to generate
%   list of obstacles
for i = 1:n    
    transform=[-1+2*rand() -1+2*rand() -limit+2*limit*rand();
        -1+2*rand() -1+2*rand() -limit+2*limit*rand()
        0 0 1];
    new=transform*[x y ones(4,1)].';
    xv(:,i)=new(1,:).';
    yv(:,i)=new(2,:).';
end

figure
plot(xv,yv);
% plot(xv,yv,1,0.5,'+');
axis equal

%Add boundaries
boundaryx=[-limit -limit limit limit;
            -limit limit limit -limit;
            -limit -limit limit limit;
            -limit limit limit -limit];
boundaryy=[-limit limit limit -limit;
            limit limit -limit -limit;
            -limit limit limit -limit;
            limit limit -limit -limit];
%xv(:,i+1)=[limit; limit; NaN*ones(m-2,1)]; yv(:,i+1)=[limit; -limit; NaN*ones(m-2,1)];

hold on
plot(boundaryx,boundaryy,'k');
title('Test Environment');
axis equal

% Now add boundaries to list of obstacles
xv(:,i+1:i+4)=boundaryx;
yv(:,i+1:i+4)=boundaryy;
