
clc
clear all
close all
file_path =  'E:\files\tqwt_code_paper\故障轴承数据\XJTU-SP\Data\XJTU-SY_Bearing_Datasets\35Hz12kN\Bearing1_2\';% 文件夹路径
%% 全寿命振动信号
csv_acc_path_list = dir(strcat(file_path,'*.csv'));%获取该文件夹中所有csv格式的文件
csv_order_name= sort_nat({csv_acc_path_list.name}); 
csv_acc_num = length(csv_acc_path_list);%获取文件总数量
if csv_acc_num > 0 %有满足条件的文件
        for j = 1:csv_acc_num %逐一读取文件
            csv_acc_name = csv_order_name(j);% 文件名
            csv_acc =  csvread(strcat(file_path,csv_acc_name{1,1}),1,0);
            csv_acc_data(:,:,j)=csv_acc;
            fprintf('%d %d %s\n',csv_acc_num,j,strcat(file_path,csv_acc_name{1,1}));% 显示正在处理的文件名
        end
end
% 合并矩阵 时间*通道
channel=2;   %信号的通道数
csv_acc_data_change=permute(csv_acc_data,[2 1 3]);
csv_acc_data=reshape(csv_acc_data_change,channel,prod(size(csv_acc_data))/channel)';

%% 全寿命振动信号的时域图
clearvars -except csv_acc_data 
figure;plot(csv_acc_data(:,1));title('水平振动信号');set(gca,'YLim',[-30 30]);
figure;plot(csv_acc_data(:,2));title('竖直振动信号');set(gca,'YLim',[-30 30]);
