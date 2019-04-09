clear all;clc;   
hold on  
% set(0,'defaultfigurecolor','w');%设置背景颜色为白色   
axis([0,4,0,1]); 
grid off   
for a=0:0.001:4; 
    x=[0.1234];  
    u=a;  
    for n=2:150      
        x(n)=u* x(n-1)*(1-x(n-1));          
    end
    for n=100:150              
        plot(a,x(n),'k','markersize',3);                
    end
end
title('\fontsize{10}Logistic映射迭代图');          
xlabel('\fontsize{10}分支参数u'),
ylabel('\fontsize{10}输出序列分布x(n)');          
set(gca,'FontSize',10); % 设置文字大小，同时影响坐标轴标注、图例、标题等。