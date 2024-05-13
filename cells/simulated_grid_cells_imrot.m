close all
clear all
clc
% simulate 10 x 10 box
global gridscale 

gridscale = 8;

x = linspace(0, 19.9, 200) ;
y = x ; 
[X,Y] = meshgrid(x,y) ; 

r0 = pi/2;
shiftx_0_1 = 1;
shifty_0_1 = 0;

r1 = pi/6;
shiftx_0_2 = 1;
shifty_0_2 = 1;

shift_xy = zeros(2, 1);

n_cells = 40*40;
shift_vec = linspace(0, 2*pi- 2*pi/sqrt(n_cells), sqrt(n_cells)) ;
[shift_vecX,shift_vecY] = meshgrid(shift_vec,shift_vec) ; 

grid_cells_env1 = zeros(n_cells, 100, 100);
grid_cells_env2 = zeros(n_cells, 100, 100);

n = 1;
for s1 = 1:sqrt(n_cells)
    for s2 = 1:sqrt(n_cells)
        shift_xy(1) = shift_vecX(s1, s2);
        shift_xy(2) = shift_vecY(s1, s2);
        grid_cell = grid_firing_rate_map(X, Y, shift_xy, r0, shiftx_0_1, shifty_0_1);
        grid_cell = grid_cell(51:150, 51:150);
        grid_cells_env1(n, :, :) = grid_cell();

        % if n<=16
        %     figure(1)
        %     subplot(4, 4, n)
        %     imagesc(grid_cell)
        % end

        grid_cell = grid_firing_rate_map(X, Y, shift_xy, r0, shiftx_0_2, shifty_0_2);
        grid_cell = imrotate(grid_cell,r1*180/pi,'bilinear','crop');
        grid_cell = grid_cell(51:150, 51:150);
        grid_cells_env2(n, :, :) = grid_cell;
        % 
        % if n<=16
        %     figure(2)
        %     subplot(4, 4, n)
        %     imagesc(grid_cell)
        % end
        n = n + 1;
    end
end

% c11 = corrcoef(reshape(grid_cells_env1, n_cells, 10000)');
% figure(10)
% imagesc(c11)
% c22 = corrcoef(reshape(grid_cells_env2, n_cells, 10000)');
% figure(20)
% imagesc(c22)
% 
% utri = triu(ones(n_cells,n_cells), 1);
% utri_c11 = triu(c11, 1);
% utri_c11 = utri_c11(utri>0);
% utri_c22 = triu(c22, 1);
% utri_c22 = utri_c22(utri>0);
% figure(30)
% plot(utri_c11, utri_c22, '*')

simulated_grid_cells = zeros(2, n_cells, 10000);
simulated_grid_cells(1, :, :) = reshape(grid_cells_env1, n_cells, 10000); 
simulated_grid_cells(2, :, :) = reshape(grid_cells_env2,  n_cells, 10000);
% 
save(['simulated_grid_cells_',num2str(n_cells),'_gridscale',num2str(gridscale),'.mat'], "simulated_grid_cells")

function Zr = rotate(theta, X, Y)
global gridscale 
    R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
    XY = [X(:) Y(:)];     % Create Matrix Of Vectors
    rotXY = XY*R'; %MULTIPLY VECTORS BY THE ROT MATRIX 
    Xr = reshape(rotXY(:,1), size(X,1), []);
    Yr = reshape(rotXY(:,2), size(Y,1), []);
    Zr = cos((2*pi/gridscale) * (Xr + Yr)) ;    
end

function grid_cell = grid_firing_rate_map(X, Y, shift_xy, r0, shiftx_0, shifty_0)
    X = X + shift_xy(1) + shiftx_0;
    Y = Y + shift_xy(2) + shifty_0;

    Z = rotate(r0, X, Y);
    theta = pi/3 + r0;
    Zr1 = rotate(theta, X, Y);
    theta = 2*pi/3 + r0;
    Zr2 = rotate(theta, X, Y);

    grid_cell = (2/3)*0.3*( Z + Zr1 + Zr2) + 0.5;
   % grid_cell(grid_cell<thresh) = 0;
end