
clc
clear all
close all
file_path =  'E:\files\tqwt_code_paper\�����������\XJTU-SP\Data\XJTU-SY_Bearing_Datasets\35Hz12kN\Bearing1_2\';% �ļ���·��
%% ȫ�������ź�
csv_acc_path_list = dir(strcat(file_path,'*.csv'));%��ȡ���ļ���������csv��ʽ���ļ�
csv_order_name= sort_nat({csv_acc_path_list.name}); 
csv_acc_num = length(csv_acc_path_list);%��ȡ�ļ�������
if csv_acc_num > 0 %�������������ļ�
        for j = 1:csv_acc_num %��һ��ȡ�ļ�
            csv_acc_name = csv_order_name(j);% �ļ���
            csv_acc =  csvread(strcat(file_path,csv_acc_name{1,1}),1,0);
            csv_acc_data(:,:,j)=csv_acc;
            fprintf('%d %d %s\n',csv_acc_num,j,strcat(file_path,csv_acc_name{1,1}));% ��ʾ���ڴ�����ļ���
        end
end
% �ϲ����� ʱ��*ͨ��
channel=2;   %�źŵ�ͨ����
csv_acc_data_change=permute(csv_acc_data,[2 1 3]);
csv_acc_data=reshape(csv_acc_data_change,channel,prod(size(csv_acc_data))/channel)';

%% ȫ�������źŵ�ʱ��ͼ
clearvars -except csv_acc_data 
figure;plot(csv_acc_data(:,1));title('ˮƽ���ź�');set(gca,'YLim',[-30 30]);
figure;plot(csv_acc_data(:,2));title('��ֱ���ź�');set(gca,'YLim',[-30 30]);