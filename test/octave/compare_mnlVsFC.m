% compare monolayer and FC

clear
clf

nb = 6;
%monolayer
pathmnl = ['../srim/mnl_3d/b' num2str(nb) '/'];
S1 = srim2mat(pathmnl);

% FC mode
path = ['../srim/b' num2str(nb) '/'];
S = srim2mat(path);

# Size of the cell
cellsize_all = [12 6 12 6 3 500 100];
cellsize = cellsize_all(nb) ;

x = S.x/10;


figure 1
subplot(1,2,1)
surf(x,x,S.Ri3.*cellsize*10e-8)
xlabel('y (nm)');
ylabel('x (nm)');
colormap('jet')
view(40,30)
title('FC mode')
subplot(1,2,2)
surf(x,x,S1.Ri3.*cellsize*10e-8)
xlabel('y (nm)');
ylabel('x (nm)');
colormap('jet')
view(40,30)
title('monolayer mode')

figure 2
clf
subplot(1,2,1)
contour(x-25000,x-35500,S.Ri3.*cellsize*10e-8,[0:0.002:0.01 0.015 0.02 0.025])
colorbar
xlabel('y (nm)');
ylabel('x (nm)');
xlim([-1 1]*7500)
ylim([-1 1]*7500)
caxis([0 0.08])
axis equal
colormap('jet')
title('FC mode')
grid
subplot(1,2,2)
contour(x-25000,x-35500,S1.Ri3.*cellsize*10e-8,[0:0.0025:0.01 0.02 0.04 0.06])
colorbar
xlabel('y (nm)');
ylabel('x (nm)');
xlim([-1 1]*7500)
ylim([-1 1]*7500)
caxis([0 0.08])
axis equal
colormap('jet')
title('monolayer mode')
grid

