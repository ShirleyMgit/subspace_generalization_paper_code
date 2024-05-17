close all
clear all
clc
% simulate 10 x 10 box
global gridscale thresh  

thresh = 0.25;
gridscale = 2;
% the actual grid scale is: gridscale*2/sqrt(3) - this is the reciprocal
% lattice

x = linspace(0, 9.9, 100) ;
y = x ; 
[X,Y] = meshgrid(x,y) ; 

r0 = -1*pi/4;
shiftx_0_1 = 0;
shifty_0_1 = 0;

shift_xy = zeros(2, 1);

r_xy = r0 + 1*pi/4;%
dx = 0.01;
max_x = 2/sqrt(3);%because of the gridscale/ actual grid scale ratio (reciprocal lattice)
max_y = max_x*sqrt(3) / 2;
dy = dx* sqrt(3) / 2;
shift_vec_x1 = linspace(0, max_x -dx , ceil(max_x/dx)) ;
shift_vec_y1 = linspace(0, max_y -dy, ceil(max_y/dy)) ;
[shift_vecX,shift_vecY] = meshgrid(shift_vec_x1,shift_vec_y1); 

cellx = length(shift_vec_x1);
celly = length(shift_vec_y1);

for y = 2:celly
    shift_vecX(y,:) = shift_vecX(y,:) + dx*(y-1)/2;
end

figure(1000)
plot(shift_vecX(:),shift_vecY(:), 'o')

R = [cos(r_xy) -sin(r_xy); sin(r_xy) cos(r_xy)];
XYshift = [shift_vecX(:) shift_vecY(:)];     
rotXYshift = XYshift*R'; 

figure(1001)
plot(rotXYshift(:, 1), rotXYshift(:, 2), 'o')

n_cells = cellx * celly;

grid_cells_env1 = zeros(n_cells, 100, 100);

n = 1;
for s1 = 1:n_cells
        shift_xy(1) = rotXYshift(s1, 1);
        shift_xy(2) = rotXYshift(s1, 2);
        grid_cell = grid_firing_rate_map(X, Y, shift_xy, r0, shiftx_0_1, shifty_0_1);
        grid_cells_env1(n, :, :) = grid_cell;

        if n<=25
            figure(1)
            subplot(5, 5, n)
            imagesc(grid_cell)
           % colorbar
            title(n)
        end

        n = n + 1;

end

mean_all1 = squeeze(mean(grid_cells_env1));

figure(100)
subplot(2,1,1)
imagesc(squeeze(grid_cells_env1(1,:,:)))
title('cell env 1')
colorbar
subplot(2,1,2)
imagesc(mean_all1)
title('mean env 1')
colorbar

c11 = corrcoef(reshape(grid_cells_env1, n_cells, 10000)');
figure(10)
imagesc(c11)
colorbar

figure(10000)
subplot(2,2,1)
imagesc(X(:), Y(:), squeeze(grid_cells_env1(1,:,:)))
colorbar
title('1')
subplot(2,2,3)
imagesc(X(:), Y(:), squeeze(grid_cells_env1(5,:,:)))
colorbar
title('5')
subplot(2,2,2)
imagesc(X(:), Y(:), squeeze(grid_cells_env1(120,:,:)))
colorbar
title('120')
subplot(2,2,4)
imagesc(X(:), Y(:), squeeze(grid_cells_env1(1,:,:)) + squeeze(grid_cells_env1(120,:,:)))
colorbar
title('1 + 120')


function Zr = rotate(theta, X, Y)
global gridscale 
    R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
    XY = [X(:) Y(:)];     
    rotXY = XY*R';
    Xr = reshape(rotXY(:,1), size(X,1), []);
    Yr = reshape(rotXY(:,2), size(Y,1), []);
    Zr = cos((2*pi/gridscale) * ((Xr + Yr)/sqrt(2))) ;    
end

function grid_cell = grid_firing_rate_map(X, Y, shift_xy, r0, shiftx_0, shifty_0)
global gridscale thresh
    X = X + (shift_xy(1) + shiftx_0)*(gridscale);
    Y = Y + (shift_xy(2) + shifty_0)*(gridscale);

    Z = rotate(r0, X, Y);

    theta = pi/3 + r0;
    Zr1 = rotate(theta, X, Y);

    theta = 2*pi/3 + r0;
    Zr2 = rotate(theta, X, Y);

    grid_cell = (1/3)*(Z + Zr1 + Zr2);% + 0.5;%0.5*( Z + Zr1); %
    grid_cell(grid_cell<thresh) = 0;
end