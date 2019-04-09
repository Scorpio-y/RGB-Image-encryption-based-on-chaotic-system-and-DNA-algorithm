clc
clear

[T,Y]=ode45(@chao_SimpleLorenz,0:0.01:500,[0.1;0.2;0.3;0.4]);
maxX = ceil(max(Y(:,1)));
minX = floor(min(Y(:,1)));
maxY = ceil(max(Y(:,2)));
minY = floor(min(Y(:,2)));
maxZ = ceil(max(Y(:,3)));
minZ = floor(min(Y(:,3)));
maxH = ceil(max(Y(:,4)));
minH = floor(min(Y(:,4)));

figure
plot3(Y(10001:end,1),Y(10001:end,2),Y(10001:end,3))
xlim([minX maxX])
ylim([minY maxY])
zlim([minZ maxZ])
xlabel('\itx')
ylabel('\ity')
zlabel('\itz')

figure
plot3(Y(10001:end,1),Y(10001:end,2),Y(10001:end,4))
xlim([minX maxX])
ylim([minY maxY])
zlim([minH maxH])
xlabel('\itx')
ylabel('\ity')
zlabel('\ith')

figure
plot3(Y(10001:end,1),Y(10001:end,3),Y(10001:end,4))
xlim([minX maxX])
ylim([minZ maxZ])
zlim([minH maxH])
xlabel('\itx')
ylabel('\itz')
zlabel('\ith')

figure
plot3(Y(10001:end,2),Y(10001:end,3),Y(10001:end,4))
xlim([minY maxY])
ylim([minZ maxZ])
zlim([minH maxH])
xlabel('\ity')
ylabel('\itz')
zlabel('\ith')