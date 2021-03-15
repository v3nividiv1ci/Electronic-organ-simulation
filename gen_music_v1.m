fs=5000; %数字信号的采样率
%s0=[0;0;0;0]; %输出数据变量初始化
s0 = [];
%s0 = 0;
D=4; %控制音乐时长系数
a=4; %控制音乐类型（包络衰减）系数
Z=2^(1/12); %音阶的半音频率系数
tone=220*Z^10; do=tone; re=tone*Z^2; mi=tone*Z^4; fa=tone*Z^5;
so=tone*Z^7; la=tone*Z^9; xi=tone*Z^11;lowre=220;lowmi=220*Z^2;
lowfa=220*Z^3; lowso=220*Z^5;lowla=220*Z^7;lowxi=220*Z^9;none=100000;
shenfa=220*Z^4;shenso=220*Z^6;
% 示例的唱名和音长
M=[lowla lowxi do lowxi do mi lowxi none lowmi lowla lowso lowla do lowso ....
    none lowmi lowmi lowfa lowmi lowfa do lowmi none do do lowxi shenfa lowfa ...
    lowxi lowxi none lowla lowxi do lowxi do mi lowxi none lowmi lowla lowso...
    lowla do lowso none lowre lowmi lowfa do lowxi do re re mi do do lowxi...
    lowla lowla lowxi shenso lowla none do re mi re mi so re none do do lowxi...
    do mi mi none lowla lowxi do lowxi do re do lowso lowso fa mi re do mi ...
    mi mi la la so so mi re do do re do re so mi mi la la so so mi re do do ...
    re do re lowxi lowla lowla];
L=[1 1 3 1 2 2 4 2 2 3 1 2 2 4 2 1 1 3 1 2 2 4 2 1 1 3 1 2 2 4 2 1 1 3 1 2 ...
    2 4 2 2 3 1 2 2 4 2 1 1 2 1 2 2 1 3 1 4 1 1 1 1 2 2 4 2 1 1 3 1 2 2 6 2 ...
    2 1 1 2 2 6 2 1 1 2 1 1 2 3 1 4 2 2 2 1 1 6 2 3 1 3 1 1 1 4 2 3 1 2 2 6 ... 
    2 3 1 3 1 1 1 4 2 3 1 2 1.5 0.5 4 ]/D;

for i=1:4
    s0i = 0;
    if(i==2 || i==4)
        A =[1 0.2 0.3 0.25 0.25 0.2]; % 各次谐波的幅度强度分布
    else
        A=[1 0 0 0 0 0]; % 无谐波
    end
    t=0:1/fs:2; 
    if(i==1 || i==2)
        B=exp(0*t); % 无包络
    else
        B=exp(-a*t); % 指数函数形式包络
        % 一次函数形式包络：B=a*(t-2*pi);
    end
    C=[0.2 0.6 1.5 2.2 0.9 1.3]; % 各次谐波的初始相位
    for k=1:length(M) % 逐个音符生成数字音乐的数据
        s=0;
        t=0: 1/fs: L(k); % 每个音符时间长度数据
        for j=1:length(A) % 逐次谐波生成数字音乐的数据
            s=s+A(j)*cos(2*pi*M(k)*j*t+C(j)); % 逐次谐波正弦信号合成
        end
        %s = cos(2*pi*M(k)*t);
        s=s .* B(1:length(s)); % 逐个音符的音乐形式（包络）修正
        %s(i,:)=zeros(1,length(s));
        s0i=[s0i,s(2:end)]; % 逐个音符的数据长度修正
    end
    s0 = [s0; s0i];
end


% 分析包络
% --------------------------------------------------------------------------
s04=s0(4,:);
format = abs(s04(1:floor(fs*1/D))); % 取绝对值, 丰富正半轴的包络
Bm = zeros(1, length(format)); % 初始化存放包络的向量
% w = floor(1 / 3 * length(format) * M(1)/fs); %这里的1/3没有特别的含义，只是测出来结果比较好看
% 约等于 27 
w = 30;
b0 = 1;
while b0+w < length(format) %平移窗口获得极大值和位置
    [bm,bx]=max(format(b0:b0+w));
    Bm(b0+bx-1)=bm;
    b0=b0+floor(w/2);
end
% Bm=Bm/max(Bm); %最大值归一化
Bm(Bm==0)=[]; %去空隙抽取包络
t=0 :length(Bm)/fs*D: length(Bm); %包络插值为1秒数据
Binter=interp1(Bm,t,'spline');

% 取消上面注释即可画图
% figure(4);scatter((1: length(Bm)), Bm); title('极值点');
% --------------------------------------------------------------------------
% 获得谐波
s0fft = abs(fft(s0(4,:)));
xf = s0fft(2:length(s0fft)); % 选择分析数据, 去除直流分量
Am = linspace(0, 0, length(xf)); % 频率的幅度
fk = linspace(0, 0, length(xf)); % 频率的位置
[Am(1),f0]=max(xf); %确定基频f0幅度和位置
fk(1)=f0;
k=2;
w=5;
while f0*k+w < length(xf) %开窗确定谐频幅度和位置
    [Am(k),maxf]=max(xf(f0*k-w:f0*k+w));
    fk(k)=f0*k-w+maxf-1;
    k=k+1; % 检索每个基频的倍数附件的值
end
Am=Am/Am(1); %归一化数据
figure;
    subplot(2, 2, 1); % 2x2 图表中的第一个图
        plot(format, 'green'); title('第一个音节的时域图');
    subplot(2, 2, 3); % 2x2 图表中的第三个图
        plot(format, 'green'); hold on;
        plot(Binter, 'red'); hold off;
        title('包络线');
    subplot(2, 2, 2); % 2x2 图表中的第二个图         
        plot(s0fft);
        title('频域图');
    subplot(2, 2, 4);% 2x2 图表中的第四个图   
        stem(fk, Am);
        title('谐波图');

Y = s0(4,:);
sound(Y,fs)   %播放语音
figure(2)
plot(Y)          %波形图
title('原始语音信号')
% grid on;

title0=['无谐波无包络语谱图'; '有谐波无包络语谱图' ;'无谐波有包络语谱图'; '有谐波有包络语谱图'];
figure(3);
for i=1:4
    subplot(2,2,i);
    Y=s0(i,:);
    spectrogram(Y,256,128,256,41000,'yaxis');
    xlabel('时间(s)')
    ylabel('频率(Hz)')
    title(title0(i,:))
end
audiowrite('try0.wav',s0(4,:),fs);