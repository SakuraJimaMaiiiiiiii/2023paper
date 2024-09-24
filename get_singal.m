%% 读取文件
flstr= dir('E:\files\tqwt_code_data\故障轴承数据\XJTU-SP\Data\XJTU-SY_Bearing_Datasets\35Hz12kN\Bearing1_2\');                   %(file root directory),获取文件名
flname1=struct2cell(flstr); 
flname2=flname1(1,:);
flnumber=size(flname2,2); 
flname_arry=flname2(3:flnumber);
m=length(flname_arry);
for i=1:m
    flname=strcat('E:\files\tqwt_code_paper\故障轴承数据\XJTU-SP\Data\XJTU-SY_Bearing_Datasets\35Hz12kN\Bearing1_2\',num2str(i),'.csv');       %获得文件全地址
    datain = csvread(flname,1,0);                                               %datain是读入的振动数据矩阵
    data_ho(:,i)=datain(:,1);
    data_ve(:,i)=datain(:,2);
    disp(num2str(i));
end
save('XJTUbearingdata.mat','data_ho','data_ve')



%% data input and display
load ('XJTUbearingdata.mat');
[len_m,len_n]=size(data_ho);
for i=1:len_n
    rms_ho(i)=rms(data_ho(:,i));
    rms_ve(i)=rms(data_ve(:,i));
end
figure,plot(rms_ho)
xlabel('时间 [分钟]'),ylabel('RMS')
figure,plot(rms_ve,'b')
xlabel('time [min]'),ylabel('Amplitude')