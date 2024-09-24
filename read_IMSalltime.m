
clc
clear all
close all
delimiter = ',';
file_path =  'E:\files\tqwt_code_paper\故障轴承数据\NASA\2nd_test\';% 文件夹路径
txt_path_list = dir(strcat(file_path,'*.txt'));%获取该文件夹中所有txt格式的文件
txt_num = length(txt_path_list);%获取文件总数量
channel=4;        %信号的通道数,IMS一共3个实验，第一个是8通道，其他都是4通道
text_data=zeros(20480,channel,txt_num);
if txt_num > 0 %有满足条件的文件
        for j = 1:txt_num %逐一读取文件
            text_name = txt_path_list(j).name;% 文件名
            text_fileID = fopen(strcat(file_path,text_name),'r');
            text = cell2mat(textscan(text_fileID,'%f%f%f%f', 'Delimiter', delimiter));
            % 注意上面那个‘%f%f%f%f’也要随通道数改变，8个就是‘%f%f%f%f%f%f%f%f’，4个就是'%f%f%f%f'
            fclose(text_fileID );
            text_data(:,:,j)=text;
            fprintf('%d %d %s\n',i,j,strcat(file_path,text_name));% 显示正在处理的文件名
        end
end
% 合并矩阵 时间*通道
text_data_change=permute(text_data,[2 1 3]);
text_data=reshape(text_data_change,channel,prod(size(text_data))/channel)';
clearvars -except text_data

figure,plot(text_data,'color', [0  114  189]/255)
xlim([0,2*10^7])
xlabel('时间[10^4s]')
ylabel('幅值')
get(gca)
set(gca, 'XTickLabel', 0 :  984)



