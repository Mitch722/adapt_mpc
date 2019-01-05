
% gp_matern_1.mat
% matern52_1.mat
% gp5.mat

load('matern52_1.mat')

[carbon, idx] = sort(y_pred);

xxx = x_pred(idx, 1);
yyy = x_pred(idx, 2);

xx2 = x_pred(idx, 3);
yy2 = x_pred(idx, 4); 

zzz = y_pred(idx);

cmp=jet(numel(y_pred));

figure;
scatter3(xxx, yyy, zzz, [], zzz, 'filled'); 
grid on;

colorbar;

figure;
scatter3(xx2, yy2, zzz, [], zzz, 'filled'); 
grid on;

colorbar;


tri = delaunay(xxx,yyy);
figure;
plot3(xxx,yyy,zzz,'.-')
plot(xxx,yyy,'.')
%
% How many triangles are there?
[r,c] = size(tri);
disp(r)
% Plot it with TRISURF
h = trisurf(tri, xxx, yyy, zzz);
axis vis3d
% Clean it up



tri = delaunay(xx2,yy2);
figure;
plot3(xx2,yy2,zzz,'.-')
plot(xx2,yy2,'.')
%
% How many triangles are there?
[r,c] = size(tri);
disp(r)
% Plot it with TRISURF
h = trisurf(tri, xx2, yy2, zzz);
axis vis3d
% Clean it up


