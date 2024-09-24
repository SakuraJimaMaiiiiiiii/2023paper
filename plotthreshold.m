clc; 
clear;
x=-0.8:0.01:0.8;
y1=0.*(abs(x)<=0.5*0.4)+(x.*(abs(x)-0.5*0.4)/(0.4/0.5-0.4*0.5)).*(abs(x)>=0.5*0.4&abs(x)<0.4/0.5)+x.*(abs(x)>=0.4/0.5);
plot(x,y1,'r','linewidth',2)
legend('改进阈值函数')
set(gca, 'XAxisLocation', 'origin')
set(gca, 'YAxisLocation', 'origin')
grid on
% set(gca, 'XAxisLocation', 'origin')
% set(gca, 'YAxisLocation', 'origin')
% set(gca,'ytick',[])  %隐去y轴坐标值
% set(gca,'xtick',[])  %隐去x轴坐标值

hold on
y2 = x.*(abs(x)>=0.4&abs(x)<=0.8)+0*(abs(x)<=0.4); %hard 
plot(x,y2,'b','linewidth',2)
% legend('硬阈值函数')
hold 
y3=(0.5*x-0.2).*(x>=0.4&x<=0.8)+0*(abs(x)<=0.4)+(0.5*x+0.2).*(x>=-0.8&x<=-0.4);
plot(x,y3,'y','linewidth',2)
% legend('软阈值函数')
legend('改进阈值函数','硬阈值函数','软阈值函数','Location','northwest')
axis tight
grid on