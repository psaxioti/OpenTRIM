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
contour(x,x,S.Ri3.*cellsize*10e-8,100)
colorbar
xlabel('y (nm)');
ylabel('x (nm)');
colormap('jet')
title('FC mode')
subplot(1,2,2)
contour(x,x,S1.Ri3.*cellsize*10e-8,100)
colorbar
xlabel('y (nm)');
ylabel('x (nm)');
colormap('jet')
title('monolayer mode')

