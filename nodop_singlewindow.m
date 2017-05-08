%% LTE上行随机接入前导检测过程
%小区半径200km，系统带宽20M，采样率30.72MHz

clear all
f=1250;                     %子载波间隔
Nzc=839;                    %ZC序列长度
Ncp=21024;                %cp长度
Nseq=24576;              %1个序列长度，在格式3中发送的为2个序列
Ngt=40960;                  %保护时间
Ncs=0;                          %循环移位值
capa=1;                         %每一个根序列可以产生的ZC序列个数
user_num=1;                  %接入数量
user_seq=randperm(capa);                %接入用户的序列
u=133;                           %ZC序列根序号
Nlocal=1024;                    %本地用于检测序列的长度
down_sample=Nseq/Nlocal;                 %降采样倍数24
snr_set=-24:-5;                          %信噪比集
samesnr=100000;

fsys=30.72*10^6;

t_delay=0.0013;                    %延时时间1.3ms
fs=f*Nlocal;                        %采样频率2.56MHz
delaymax=t_delay*fs;            %延时采样点数
du=100;                             %对应根序列的du
couter=0;

right=zeros(1,length(snr_set));
false_counter=zeros(1,length(snr_set));
miss_counter=zeros(1,length(snr_set));
%% 前导产生
ZC_root=zeros(1,Nzc);                           %根序列
for n=0:Nzc-1
    ZC_root(1,n+1)=exp(-1i*pi*u*n*(n+1)/Nzc);
end

%进行循环移位并dft
Cv=0;                               %循环移位值
ZCf_root=1/sqrt(Nzc)*fft(ZC_root);
ZC=ZC_root;
ZCf=1/sqrt(Nzc)*fft(ZC);


%% 产生本地用于检测的ZC序列
ZCf_zero=zeros(1,Nlocal);
ZCf_zero(floor((Nlocal-Nzc)/2)+1:floor((Nlocal+Nzc)/2))=ZCf_root;                  %频域补零，补零补在高频
ZCf_local=conj(ZCf_zero);

%% 在不同信噪比下发送preamble
for snrloop=1:length(snr_set)
    snrrev=snr_set(snrloop);
   
for timeloop=1:samesnr        %相同的snr多次运算得出虚警率
    
    %把序列放到频谱对应位置上
    ZCf_location=[zeros(1,12),ZCf,zeros(1,13),zeros(1,Nseq-Nzc-25)];                                     %放在频域的最前面6个RB[(1,864),(1,23712)]
    ZCf_locaf=[ZCf_location(floor(Nseq/2)+1:end),ZCf_location(1:floor(Nseq/2))];                    %为了抵消fft对频谱的搬移作用
    
    %ifft回时域
    ZC_ifft=sqrt(Nseq)*ifft(ZCf_locaf);
    
    %两个序列长度
    ZC_double=[ZC_ifft,ZC_ifft];
    
    %加cp
    ZC_cp=[ZC_double(end-Ncp+1:end),ZC_double];

    %加保护时间
    ZC_safe=[ZC_cp,zeros(1,Ngt)];
    
    %加传输延时
    delay=down_sample*randperm(delaymax);                    %每个用户随机延迟
    delay_down=delay/down_sample;
    ZC_delay=[zeros(1,delay(1)),ZC_cp,zeros(1,Ngt-delay(1))];


%% 加AWGN信道
ZC_awgn=awgn(ZC_delay,snrrev-12.6761,'measured');

% 信号功率
P_signal=mean(abs(ZC_delay).^2);
% 接收信号功率
P_noise=mean(abs(ZC_awgn).^2);

%% 接收信号

%去cp
ZCrev_cp=ZC_awgn(Ncp+1:end);

    ZCrev_det=ZCrev_cp(1:Nseq);
    ZCrev_det2=ZCrev_cp(Nseq+1:Nseq*2);             %重复序列

     %将频谱解映射，移到低频进行滤波
     k=11856;                    %将频谱从最前端搬移到中间       (24576-864)/2
    m=1:Nseq;
    ZCrev_low=ZCrev_det.*exp(1i*2*pi*k*m/Nseq);
    ZCrev_low2=ZCrev_det2.*exp(1i*2*pi*k*m/Nseq);

    %经过低通滤波器
    load lowpassfilter.mat
    ZCrev_filter=filter(lowpassfilter,1,ZCrev_low);
    ZCrev_filter2=filter(lowpassfilter,1,ZCrev_low2);
    
    %进行降采样
    ZCrev_down=downsample(ZCrev_filter,down_sample);
    ZCrev_down2=downsample(ZCrev_filter2,down_sample);
    
    %进行fft到频域
    ZCrev_fft=1/sqrt(Nlocal)*fft(ZCrev_down);                                   %序列1和序列2相当于在时域移了2048个点
    ZCrev_fftshift=[ZCrev_fft(floor(Nlocal/2)+1:end),ZCrev_fft(1:floor(Nlocal/2))];
    ZCrev_fft2=1/sqrt(Nlocal)*fft(ZCrev_down2);
    ZCrev_fftshift2=[ZCrev_fft2(floor(Nlocal/2)+1:end),ZCrev_fft2(1:floor(Nlocal/2))];
    
    %频域共轭相乘
    ZCrev_mul=ZCrev_fftshift.*ZCf_local;
    ZCrev_mul2=ZCrev_fftshift2.*ZCf_local;
    
    ZCrev_time=sqrt(Nlocal)*ifft(ZCrev_mul);
    ZCrev_real=abs(real(ZCrev_time));
    ZCrev_imag=abs(imag(ZCrev_time));
    ZCrev_pow=(max(ZCrev_real,ZCrev_imag)+1/2*min(ZCrev_real,ZCrev_imag)).^2;
    
    ZCrev_time2=sqrt(Nlocal)*ifft(ZCrev_mul2);
    ZCrev_real2=abs(real(ZCrev_time2));
    ZCrev_imag2=abs(imag(ZCrev_time2));
    ZCrev_pow2=(max(ZCrev_real2,ZCrev_imag2)+1/2*min(ZCrev_real2,ZCrev_imag2)).^2;
    
    [maxval,location]=max(ZCrev_pow+ZCrev_pow2);
    
      %% 设置门限
menxian=20;       %ZC序列的门限值
    if location==1
        reallocation=Ncp/24+location-2;
    else
        if maxval>menxian
            reallocation=location-2;               %滤波器群时延和序列第一位是从1开始而不是0引起的
        end
    end
    
    
    if maxval<menxian
        miss_counter(snrloop)=miss_counter(snrloop)+1;
        reallocation=0;
    end


if reallocation==delay_down(1)
    right(snrloop)=right(snrloop)+1;
end

end

rightdete(snrloop)=right(snrloop)/samesnr;
end

%% 画图
plot(snr_set,rightdete,'-ro')
axis([-24 -5 0 1]);
grid on;
title('正确检测率(AWGN信道)');
xlabel('SNR(dB)');
ylabel('正确检测率');
    






