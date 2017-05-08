%% LTE�����������ǰ��������
%С���뾶200km��ϵͳ����20M��������30.72MHz

clear all
f=1250;                     %���ز����
Nzc=839;                    %ZC���г���
Ncp=21024;                %cp����
Nseq=24576;              %1�����г��ȣ��ڸ�ʽ3�з��͵�Ϊ2������
Ngt=40960;                  %����ʱ��
Ncs=0;                          %ѭ����λֵ
capa=1;                         %ÿһ�������п��Բ�����ZC���и���
user_num=1;                  %��������
user_seq=randperm(capa);                %�����û�������
u=133;                           %ZC���и����
Nlocal=1024;                    %�������ڼ�����еĳ���
down_sample=Nseq/Nlocal;                 %����������24
snr_set=-24:-5;                          %����ȼ�
samesnr=100000;

fsys=30.72*10^6;

t_delay=0.0013;                    %��ʱʱ��1.3ms
fs=f*Nlocal;                        %����Ƶ��2.56MHz
delaymax=t_delay*fs;            %��ʱ��������
du=100;                             %��Ӧ�����е�du
couter=0;

right=zeros(1,length(snr_set));
false_counter=zeros(1,length(snr_set));
miss_counter=zeros(1,length(snr_set));
%% ǰ������
ZC_root=zeros(1,Nzc);                           %������
for n=0:Nzc-1
    ZC_root(1,n+1)=exp(-1i*pi*u*n*(n+1)/Nzc);
end

%����ѭ����λ��dft
Cv=0;                               %ѭ����λֵ
ZCf_root=1/sqrt(Nzc)*fft(ZC_root);
ZC=ZC_root;
ZCf=1/sqrt(Nzc)*fft(ZC);


%% �����������ڼ���ZC����
ZCf_zero=zeros(1,Nlocal);
ZCf_zero(floor((Nlocal-Nzc)/2)+1:floor((Nlocal+Nzc)/2))=ZCf_root;                  %Ƶ���㣬���㲹�ڸ�Ƶ
ZCf_local=conj(ZCf_zero);

%% �ڲ�ͬ������·���preamble
for snrloop=1:length(snr_set)
    snrrev=snr_set(snrloop);
   
for timeloop=1:samesnr        %��ͬ��snr�������ó��龯��
    
    %�����зŵ�Ƶ�׶�Ӧλ����
    ZCf_location=[zeros(1,12),ZCf,zeros(1,13),zeros(1,Nseq-Nzc-25)];                                     %����Ƶ�����ǰ��6��RB[(1,864),(1,23712)]
    ZCf_locaf=[ZCf_location(floor(Nseq/2)+1:end),ZCf_location(1:floor(Nseq/2))];                    %Ϊ�˵���fft��Ƶ�׵İ�������
    
    %ifft��ʱ��
    ZC_ifft=sqrt(Nseq)*ifft(ZCf_locaf);
    
    %�������г���
    ZC_double=[ZC_ifft,ZC_ifft];
    
    %��cp
    ZC_cp=[ZC_double(end-Ncp+1:end),ZC_double];

    %�ӱ���ʱ��
    ZC_safe=[ZC_cp,zeros(1,Ngt)];
    
    %�Ӵ�����ʱ
    delay=down_sample*randperm(delaymax);                    %ÿ���û�����ӳ�
    delay_down=delay/down_sample;
    ZC_delay=[zeros(1,delay(1)),ZC_cp,zeros(1,Ngt-delay(1))];


%% ��AWGN�ŵ�
ZC_awgn=awgn(ZC_delay,snrrev-12.6761,'measured');

% �źŹ���
P_signal=mean(abs(ZC_delay).^2);
% �����źŹ���
P_noise=mean(abs(ZC_awgn).^2);

%% �����ź�

%ȥcp
ZCrev_cp=ZC_awgn(Ncp+1:end);

    ZCrev_det=ZCrev_cp(1:Nseq);
    ZCrev_det2=ZCrev_cp(Nseq+1:Nseq*2);             %�ظ�����

     %��Ƶ�׽�ӳ�䣬�Ƶ���Ƶ�����˲�
     k=11856;                    %��Ƶ�״���ǰ�˰��Ƶ��м�       (24576-864)/2
    m=1:Nseq;
    ZCrev_low=ZCrev_det.*exp(1i*2*pi*k*m/Nseq);
    ZCrev_low2=ZCrev_det2.*exp(1i*2*pi*k*m/Nseq);

    %������ͨ�˲���
    load lowpassfilter.mat
    ZCrev_filter=filter(lowpassfilter,1,ZCrev_low);
    ZCrev_filter2=filter(lowpassfilter,1,ZCrev_low2);
    
    %���н�����
    ZCrev_down=downsample(ZCrev_filter,down_sample);
    ZCrev_down2=downsample(ZCrev_filter2,down_sample);
    
    %����fft��Ƶ��
    ZCrev_fft=1/sqrt(Nlocal)*fft(ZCrev_down);                                   %����1������2�൱����ʱ������2048����
    ZCrev_fftshift=[ZCrev_fft(floor(Nlocal/2)+1:end),ZCrev_fft(1:floor(Nlocal/2))];
    ZCrev_fft2=1/sqrt(Nlocal)*fft(ZCrev_down2);
    ZCrev_fftshift2=[ZCrev_fft2(floor(Nlocal/2)+1:end),ZCrev_fft2(1:floor(Nlocal/2))];
    
    %Ƶ�������
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
    
      %% ��������
menxian=20;       %ZC���е�����ֵ
    if location==1
        reallocation=Ncp/24+location-2;
    else
        if maxval>menxian
            reallocation=location-2;               %�˲���Ⱥʱ�Ӻ����е�һλ�Ǵ�1��ʼ������0�����
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

%% ��ͼ
plot(snr_set,rightdete,'-ro')
axis([-24 -5 0 1]);
grid on;
title('��ȷ�����(AWGN�ŵ�)');
xlabel('SNR(dB)');
ylabel('��ȷ�����');
    






