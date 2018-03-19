%% 基于混沌系统与DNA编码的彩色数字图像加密系统
clear;clc;
I=imread('../原始、加密、解密图片/lena.png','png');         %读取图像信息
I1=I(:,:,1);        %R
I2=I(:,:,2);        %G
I3=I(:,:,3);        %B
figure;imshow(I);title('原始图片');
figure;imhist(I1);title('原始图片R通道直方图');
figure;imhist(I2);title('原始图片G通道直方图');
figure;imhist(I3);title('原始图片B通道直方图');
% axis([0 255 0 4000]);
[M,N]=size(I1);                      %将图像的行列赋值给M,N
t=4;    %分块大小

%% 原始图片R,G,B通道信息熵
%R通道
T1_R=imhist(I1);   %统计图像R通道灰度值从0~255的分布情况，存至T1
S1_R=sum(T1_R);     %计算整幅图像R通道的灰度值
xxs1_R=0;           %原始图片R通道相关性
%G通道
T1_G=imhist(I2);
S1_G=sum(T1_G);
xxs1_G=0;
%B通道
T1_B=imhist(I3);
S1_B=sum(T1_B);
xxs1_B=0;

for i=1:256
    pp1_R=T1_R(i)/S1_R;   %每个灰度值占比，即每个灰度值的概率
    pp1_G=T1_G(i)/S1_G;
    pp1_B=T1_B(i)/S1_B;
    if pp1_R~=0
        xxs1_R=xxs1_R-pp1_R*log2(pp1_R);
    end
    if pp1_G~=0
        xxs1_G=xxs1_G-pp1_G*log2(pp1_G);
    end
    if pp1_B~=0
        xxs1_B=xxs1_B-pp1_B*log2(pp1_B);
    end
end

%% 原始图像相邻像素相关性分析
%{
先随机在0~M-1行和0~N-1列选中1000个像素点，
计算水平相关性时，选择每个点的相邻的右边的点；
计算垂直相关性时，选择每个点的相邻的下方的点；
计算对角线相关性时，选择每个点的相邻的右下方的点。
%}
NN=1000;    %随机取1000对像素点
x1=ceil(rand(1,NN)*(M-1));      %生成1000个1~M-1的随机整数作为行
y1=ceil(rand(1,NN)*(N-1));      %生成1000个1~N-1的随机整数作为列
%预分配内存
XX_R_SP=zeros(1,1000);YY_R_SP=zeros(1,1000);        %水平
XX_G_SP=zeros(1,1000);YY_G_SP=zeros(1,1000);
XX_B_SP=zeros(1,1000);YY_B_SP=zeros(1,1000);
XX_R_CZ=zeros(1,1000);YY_R_CZ=zeros(1,1000);        %垂直
XX_G_CZ=zeros(1,1000);YY_G_CZ=zeros(1,1000);
XX_B_CZ=zeros(1,1000);YY_B_CZ=zeros(1,1000);
XX_R_DJX=zeros(1,1000);YY_R_DJX=zeros(1,1000);      %对角线
XX_G_DJX=zeros(1,1000);YY_G_DJX=zeros(1,1000);
XX_B_DJX=zeros(1,1000);YY_B_DJX=zeros(1,1000);
for i=1:1000
    %水平
    XX_R_SP(i)=I1(x1(i),y1(i));
    YY_R_SP(i)=I1(x1(i)+1,y1(i));
    XX_G_SP(i)=I2(x1(i),y1(i));
    YY_G_SP(i)=I2(x1(i)+1,y1(i));
    XX_B_SP(i)=I3(x1(i),y1(i));
    YY_B_SP(i)=I3(x1(i)+1,y1(i));
    %垂直
    XX_R_CZ(i)=I1(x1(i),y1(i));
    YY_R_CZ(i)=I1(x1(i),y1(i)+1);
    XX_G_CZ(i)=I2(x1(i),y1(i));
    YY_G_CZ(i)=I2(x1(i),y1(i)+1);
    XX_B_CZ(i)=I3(x1(i),y1(i));
    YY_B_CZ(i)=I3(x1(i),y1(i)+1);
    %对角线
    XX_R_DJX(i)=I1(x1(i),y1(i));
    YY_R_DJX(i)=I1(x1(i)+1,y1(i)+1);
    XX_G_DJX(i)=I2(x1(i),y1(i));
    YY_G_DJX(i)=I2(x1(i)+1,y1(i)+1);
    XX_B_DJX(i)=I3(x1(i),y1(i));
    YY_B_DJX(i)=I3(x1(i)+1,y1(i)+1);
end
%水平
figure;scatter(XX_R_SP,YY_R_SP,18,'filled');xlabel('R通道随机点像素灰度值');ylabel('与该点相邻水平方向像素灰度值');title('原始图像R通道水平相关性曲线');
figure;scatter(XX_G_SP,YY_G_SP,18,'filled');xlabel('G通道随机点像素灰度值');ylabel('与该点相邻水平方向像素灰度值');title('原始图像G通道水平相关性曲线');
figure;scatter(XX_B_SP,YY_B_SP,18,'filled');xlabel('B通道随机点像素灰度值');ylabel('与该点相邻水平方向像素灰度值');title('原始图像B通道水平相关性曲线');
%垂直
figure;scatter(XX_R_CZ,YY_R_CZ,18,'filled');xlabel('R通道随机点像素灰度值');ylabel('与该点相邻垂直方向像素灰度值');title('原始图像R通道垂直相关性曲线');
figure;scatter(XX_G_CZ,YY_G_CZ,18,'filled');xlabel('G通道随机点像素灰度值');ylabel('与该点相邻垂直方向像素灰度值');title('原始图像G通道垂直相关性曲线');
figure;scatter(XX_B_CZ,YY_B_CZ,18,'filled');xlabel('B通道随机点像素灰度值');ylabel('与该点相邻垂直方向像素灰度值');title('原始图像B通道垂直相关性曲线');
%对角线
figure;scatter(XX_R_DJX,YY_R_DJX,18,'filled');xlabel('R通道随机点像素灰度值');ylabel('与该点相邻对角线方向像素灰度值');title('原始图像R通道对角线相关性曲线');
figure;scatter(XX_G_DJX,YY_G_DJX,18,'filled');xlabel('G通道随机点像素灰度值');ylabel('与该点相邻对角线方向像素灰度值');title('原始图像G通道对角线相关性曲线');
figure;scatter(XX_B_DJX,YY_B_DJX,18,'filled');xlabel('B通道随机点像素灰度值');ylabel('与该点相邻对角线方向像素灰度值');title('原始图像B通道对角线相关性曲线');
%R通道
EX1_R=0;EY1_SP_R=0;DX1_R=0;DY1_SP_R=0;COVXY1_SP_R=0;    %计算水平相关性时需要的变量
EY1_CZ_R=0;DY1_CZ_R=0;COVXY1_CZ_R=0;                %垂直
EY1_DJX_R=0;DY1_DJX_R=0;COVXY1_DJX_R=0;             %对角线
%G通道
EX1_G=0;EY1_SP_G=0;DX1_G=0;DY1_SP_G=0;COVXY1_SP_G=0;
EY1_CZ_G=0;DY1_CZ_G=0;COVXY1_CZ_G=0;
EY1_DJX_G=0;DY1_DJX_G=0;COVXY1_DJX_G=0;
%B通道
EX1_B=0;EY1_SP_B=0;DX1_B=0;DY1_SP_B=0;COVXY1_SP_B=0;
EY1_CZ_B=0;DY1_CZ_B=0;COVXY1_CZ_B=0;
EY1_DJX_B=0;DY1_DJX_B=0;COVXY1_DJX_B=0;

I1=double(I1);
I2=double(I2);
I3=double(I3);
for i=1:NN
    %第一个像素点的E，水平、垂直、对角线时计算得出的第一个像素点的E相同，统一用EX1表示
    EX1_R=EX1_R+I1(x1(i),y1(i)); 
    EX1_G=EX1_G+I2(x1(i),y1(i)); 
    EX1_B=EX1_B+I3(x1(i),y1(i)); 
    %第二个像素点的E，水平、垂直、对角线的E分别对应EY1_SP、EY1_CZ、EY1_DJX
    %R通道
    EY1_SP_R=EY1_SP_R+I1(x1(i),y1(i)+1);
    EY1_CZ_R=EY1_CZ_R+I1(x1(i)+1,y1(i));
    EY1_DJX_R=EY1_DJX_R+I1(x1(i)+1,y1(i)+1);
    %G通道
    EY1_SP_G=EY1_SP_G+I2(x1(i),y1(i)+1);
    EY1_CZ_G=EY1_CZ_G+I2(x1(i)+1,y1(i));
    EY1_DJX_G=EY1_DJX_G+I2(x1(i)+1,y1(i)+1);
    %B通道
    EY1_SP_B=EY1_SP_B+I3(x1(i),y1(i)+1);
    EY1_CZ_B=EY1_CZ_B+I3(x1(i)+1,y1(i));
    EY1_DJX_B=EY1_DJX_B+I3(x1(i)+1,y1(i)+1);
end
%统一在循环外除以像素点对数1000，可减少运算次数
% R通道
EX1_R=EX1_R/NN;
EY1_SP_R=EY1_SP_R/NN;
EY1_CZ_R=EY1_CZ_R/NN;
EY1_DJX_R=EY1_DJX_R/NN;
% G通道
EX1_G=EX1_G/NN;
EY1_SP_G=EY1_SP_G/NN;
EY1_CZ_G=EY1_CZ_G/NN;
EY1_DJX_G=EY1_DJX_G/NN;
% B通道
EX1_B=EX1_B/NN;
EY1_SP_B=EY1_SP_B/NN;
EY1_CZ_B=EY1_CZ_B/NN;
EY1_DJX_B=EY1_DJX_B/NN;
for i=1:NN
    %第一个像素点的D，水平、垂直、对角线时计算得出第一个像素点的D相同，统一用DX表示
    DX1_R=DX1_R+(I1(x1(i),y1(i))-EX1_R)^2;
    DX1_G=DX1_G+(I2(x1(i),y1(i))-EX1_G)^2;
    DX1_B=DX1_B+(I3(x1(i),y1(i))-EX1_B)^2;
    %第二个像素点的D，水平、垂直、对角线的D分别对应DY1_SP、DY1_CZ、DY1_DJX
    %R通道
    DY1_SP_R=DY1_SP_R+(I1(x1(i),y1(i)+1)-EY1_SP_R)^2;
    DY1_CZ_R=DY1_CZ_R+(I1(x1(i)+1,y1(i))-EY1_CZ_R)^2;
    DY1_DJX_R=DY1_DJX_R+(I1(x1(i)+1,y1(i)+1)-EY1_DJX_R)^2;
    %G通道
    DY1_SP_G=DY1_SP_G+(I2(x1(i),y1(i)+1)-EY1_SP_G)^2;
    DY1_CZ_G=DY1_CZ_G+(I2(x1(i)+1,y1(i))-EY1_CZ_G)^2;
    DY1_DJX_G=DY1_DJX_G+(I2(x1(i)+1,y1(i)+1)-EY1_DJX_G)^2;
    %B通道
    DY1_SP_B=DY1_SP_B+(I3(x1(i),y1(i)+1)-EY1_SP_B)^2;
    DY1_CZ_B=DY1_CZ_B+(I3(x1(i)+1,y1(i))-EY1_CZ_B)^2;
    DY1_DJX_B=DY1_DJX_B+(I3(x1(i)+1,y1(i)+1)-EY1_DJX_B)^2;
    %两个相邻像素点相关函数的计算，水平、垂直、对角线
    %R通道
    COVXY1_SP_R=COVXY1_SP_R+(I1(x1(i),y1(i))-EX1_R)*(I1(x1(i),y1(i)+1)-EY1_SP_R);
    COVXY1_CZ_R=COVXY1_CZ_R+(I1(x1(i),y1(i))-EX1_R)*(I1(x1(i)+1,y1(i))-EY1_CZ_R);
    COVXY1_DJX_R=COVXY1_DJX_R+(I1(x1(i),y1(i))-EX1_R)*(I1(x1(i)+1,y1(i)+1)-EY1_DJX_R);
    %G通道
    COVXY1_SP_G=COVXY1_SP_G+(I2(x1(i),y1(i))-EX1_G)*(I2(x1(i),y1(i)+1)-EY1_SP_G);
    COVXY1_CZ_G=COVXY1_CZ_G+(I2(x1(i),y1(i))-EX1_G)*(I2(x1(i)+1,y1(i))-EY1_CZ_G);
    COVXY1_DJX_G=COVXY1_DJX_G+(I2(x1(i),y1(i))-EX1_G)*(I2(x1(i)+1,y1(i)+1)-EY1_DJX_G);
    %B通道
    COVXY1_SP_B=COVXY1_SP_B+(I3(x1(i),y1(i))-EX1_B)*(I3(x1(i),y1(i)+1)-EY1_SP_B);
    COVXY1_CZ_B=COVXY1_CZ_B+(I3(x1(i),y1(i))-EX1_B)*(I3(x1(i)+1,y1(i))-EY1_CZ_B);
    COVXY1_DJX_B=COVXY1_DJX_B+(I3(x1(i),y1(i))-EX1_B)*(I3(x1(i)+1,y1(i)+1)-EY1_DJX_B);
end
%统一在循环外除以像素点对数1000，可减少运算次数
%R通道
DX1_R=DX1_R/NN;
DY1_SP_R=DY1_SP_R/NN;
DY1_CZ_R=DY1_CZ_R/NN;
DY1_DJX_R=DY1_DJX_R/NN;
COVXY1_SP_R=COVXY1_SP_R/NN;
COVXY1_CZ_R=COVXY1_CZ_R/NN;
COVXY1_DJX_R=COVXY1_DJX_R/NN;
%G通道
DX1_G=DX1_G/NN;
DY1_SP_G=DY1_SP_G/NN;
DY1_CZ_G=DY1_CZ_G/NN;
DY1_DJX_G=DY1_DJX_G/NN;
COVXY1_SP_G=COVXY1_SP_G/NN;
COVXY1_CZ_G=COVXY1_CZ_G/NN;
COVXY1_DJX_G=COVXY1_DJX_G/NN;
%B通道
DX1_B=DX1_B/NN;
DY1_SP_B=DY1_SP_B/NN;
DY1_CZ_B=DY1_CZ_B/NN;
DY1_DJX_B=DY1_DJX_B/NN;
COVXY1_SP_B=COVXY1_SP_B/NN;
COVXY1_CZ_B=COVXY1_CZ_B/NN;
COVXY1_DJX_B=COVXY1_DJX_B/NN;
%水平、垂直、对角线的相关性
%R通道
RXY1_SP_R=COVXY1_SP_R/sqrt(DX1_R*DY1_SP_R);
RXY1_CZ_R=COVXY1_CZ_R/sqrt(DX1_R*DY1_CZ_R);
RXY1_DJX_R=COVXY1_DJX_R/sqrt(DX1_R*DY1_DJX_R);
%G通道
RXY1_SP_G=COVXY1_SP_G/sqrt(DX1_G*DY1_SP_G);
RXY1_CZ_G=COVXY1_CZ_G/sqrt(DX1_G*DY1_CZ_G);
RXY1_DJX_G=COVXY1_DJX_G/sqrt(DX1_G*DY1_DJX_G);
%B通道
RXY1_SP_B=COVXY1_SP_B/sqrt(DX1_B*DY1_SP_B);
RXY1_CZ_B=COVXY1_CZ_B/sqrt(DX1_B*DY1_CZ_B);
RXY1_DJX_B=COVXY1_DJX_B/sqrt(DX1_B*DY1_DJX_B);

%% 1.补零
%将图像的行列数都补成可以被t整除的数，t为分块的大小。
M1=mod(M,t);    %可作为固定密钥，一遍解码时可以去除补上的0
N1=mod(N,t);    %可作为固定密钥，一遍解码时可以去除补上的0
if M1~=0
    I1(M+1:M+t-M1,:)=0;
    I2(M+1:M+t-M1,:)=0;
    I3(M+1:M+t-M1,:)=0;
end
if N1~=0
    I1(:,N+1:N+t-N1)=0;
    I2(:,N+1:N+t-N1)=0;
    I3(:,N+1:N+t-N1)=0;
end
[M,N]=size(I1);  %补零后的行数和列数
SUM=M*N;

%% 2.产生Logistic混沌序列
u=3.99;     %Logistic参数μ，自定为3.99
x0=sum(I1(:))/(255*SUM);     %计算得出Logistic初值x0
x0=floor(x0*10^4)/10^4;     %保留4位小数
p=zeros(1,SUM+1000);        %预分配内存
p(1)=x0;
for i=1:SUM+999                 %进行SUM+999次循环，共得到SUM+1000点（包括初值）
    p(i+1)=u*p(i)*(1-p(i));
end
p=p(1001:length(p));            %去除前1000点，获得更好的随机性

%% 3.将p序列变换到0~255范围内整数，转换成M*N的二维矩阵R
p=mod(ceil(p*10^3),256);
R=reshape(p,N,M)';  %转成M行N列的随机矩阵R

%% 4.求解Chen氏超混沌系统
%求四个初值X0,Y0,Z0,H0
r=(M/t)*(N/t);      %r为分块个数
%求出四个初值
X0=sum(sum(bitand(I1,3)))/(3*SUM);
Y0=sum(sum(bitand(I1,12)/4))/(3*SUM);
Z0=sum(sum(bitand(I1,48)/16))/(3*SUM);
H0=sum(sum(bitand(I1,192)/64))/(3*SUM);
%保留四位小数
X0=floor(X0*10^4)/10^4;
Y0=floor(Y0*10^4)/10^4;
Z0=floor(Z0*10^4)/10^4;
H0=floor(H0*10^4)/10^4;
%根据初值，求解Chen氏超混沌系统，得到四个混沌序列
A=chen_output(X0,Y0,Z0,H0,r);   
X=A(:,1);
X=X(1502:length(X));        %去除前1501项，获得更好的随机性（求解陈氏系统的子函数多计算了1500点）
Y=A(:,2);
Y=Y(1502:length(Y));
Z=A(:,3);
Z=Z(1502:length(Z));
H=A(:,4);
H=H(1502:length(H));

%% 5.DNA编码
%X,Y分别决定I和R的DNA编码方式，有8种，1~8
%Z决定运算方式，有3种，0~2，0表示加，1表示减，2表示异或
%H表示DNA解码方式，有8种，1~8
X=mod(floor(X*10^4),8)+1;
Y=mod(floor(Y*10^4),8)+1;
Z=mod(floor(Z*10^4),3);
H=mod(floor(H*10^4),8)+1;
e=N/t;  %e表示每一行可以分为多少块

Q2=DNA_bian(fenkuai(t,R,1),Y(1));
%R通道
Q1_R=DNA_bian(fenkuai(t,I1,1),X(1));
Q_last_R=DNA_yunsuan(Q1_R,Q2,Z(1));
Q_R(1:t,1:t)=DNA_jie(Q_last_R,H(1));
%G通道
Q1_G=DNA_bian(fenkuai(t,I2,1),X(1));
Q_last_G=DNA_yunsuan(Q1_G,Q2,Z(1));
Q_G(1:t,1:t)=DNA_jie(Q_last_G,H(1));
%B通道
Q1_B=DNA_bian(fenkuai(t,I3,1),X(1));
Q_last_B=DNA_yunsuan(Q1_B,Q2,Z(1));
Q_B(1:t,1:t)=DNA_jie(Q_last_B,H(1));
for i=2:r
    Q1_R=DNA_bian(fenkuai(t,I1,i),X(i));   %对原始图像R通道每一个分块按X对应的序号进行DNA编码
    Q1_G=DNA_bian(fenkuai(t,I2,i),X(i));
    Q1_B=DNA_bian(fenkuai(t,I3,i),X(i));
    
    Q2=DNA_bian(fenkuai(t,R,i),Y(i));   %对R的每一个分块按Y对应的序号进行DNA编码
    %R通道
    Q3_R=DNA_yunsuan(Q1_R,Q2,Z(i));         %对上面两个编码好的块按Z对应的序号进行DNA运算
    Q4_R=DNA_yunsuan(Q3_R,Q_last_R,Z(i));     %运算结果在和前一块按Z对应的序号再一次进行运算，称为扩散
    Q_last_R=Q4_R;
    %G通道
    Q3_G=DNA_yunsuan(Q1_G,Q2,Z(i));
    Q4_G=DNA_yunsuan(Q3_G,Q_last_G,Z(i));
    Q_last_G=Q4_G;
    %B通道
    Q3_B=DNA_yunsuan(Q1_B,Q2,Z(i));
    Q4_B=DNA_yunsuan(Q3_B,Q_last_B,Z(i));
    Q_last_B=Q4_B;
    
    xx=floor(i/e)+1;
    yy=mod(i,e);
    if yy==0
        xx=xx-1;
        yy=e;
    end
    Q_R((xx-1)*t+1:xx*t,(yy-1)*t+1:yy*t)=DNA_jie(Q4_R,H(i));    %将每一块合并成完整的图Q
    Q_G((xx-1)*t+1:xx*t,(yy-1)*t+1:yy*t)=DNA_jie(Q4_G,H(i));
    Q_B((xx-1)*t+1:xx*t,(yy-1)*t+1:yy*t)=DNA_jie(Q4_B,H(i));
end
Q_R=uint8(Q_R);
Q_G=uint8(Q_G);
Q_B=uint8(Q_B);

%% 抗裁剪
xx0=0.3456;
xx1=0.4532;
ppx=zeros(1,M+1000);        %预分配内存
ppy=zeros(1,N+1000); 
ppx(1)=xx0;
ppy(1)=xx1;
for i=1:M+999                 %进行SUM+999次循环，共得到SUM+1000点（包括初值）
    ppx(i+1)=u*ppx(i)*(1-ppx(i));
end
for i=1:N+999                 %进行SUM+999次循环，共得到SUM+1000点（包括初值）
    ppy(i+1)=u*ppy(i)*(1-ppy(i));
end
ppx=ppx(1001:length(ppx));            %去除前1000点，获得更好的随机性
ppy=ppy(1001:length(ppy));

[~,Ux]=sort(ppx,'descend');
[~,Uy]=sort(ppy,'descend');

for i=1:M
    temp = Q_R(i,:);
    Q_R(i,:) = Q_R(Ux(i),:);
    Q_R(Ux(i),:) = temp;
    temp = Q_G(i,:);
    Q_G(i,:) = Q_G(Ux(i),:);
    Q_G(Ux(i),:) = temp;
    temp = Q_B(i,:);
    Q_B(i,:) = Q_B(Ux(i),:);
    Q_B(Ux(i),:) = temp;
end

for i=1:N
    temp = Q_R(:,i);
    Q_R(:,i) = Q_R(:,Uy(i));
    Q_R(:,Uy(i)) = temp;
    temp = Q_G(:,i);
    Q_G(:,i) = Q_G(:,Uy(i));
    Q_G(:,Uy(i)) = temp;
    temp = Q_B(:,i);
    Q_B(:,i) = Q_B(:,Uy(i));
    Q_B(:,Uy(i)) = temp;
end

figure;imhist(Q_R);title('加密后R通道直方图');
axis([0 255 0 2000]);
figure;imhist(Q_G);title('加密后G通道直方图');
axis([0 255 0 2000]);
figure;imhist(Q_B);title('加密后B通道直方图');
axis([0 255 0 2000]);
Q_jiami(:,:,1)=Q_R;
Q_jiami(:,:,2)=Q_G;
Q_jiami(:,:,3)=Q_B;
% Q=imnoise(Q,'salt & pepper',0.1);   %加入10%的椒盐噪声
imwrite(Q_jiami,'../原始、加密、解密图片/加密后的lena.png','png');        
figure;imshow(Q_jiami);title('加密后图片');

%% 加密后信息熵
%R通道
T2_R=imhist(Q_R);
S2_R=sum(T2_R);
xxs2_R=0;
%G通道
T2_G=imhist(Q_G);
S2_G=sum(T2_G);
xxs2_G=0;
%B通道
T2_B=imhist(Q_B);
S2_B=sum(T2_B);
xxs2_B=0;
for i=1:256
    pp2_R=T2_R(i)/S2_R;
    pp2_G=T2_R(i)/S2_G;
    pp2_B=T2_R(i)/S2_B;
    if pp2_R~=0
        xxs2_R=xxs2_R-pp2_R*log2(pp2_R);
    end
    if pp2_G~=0
        xxs2_G=xxs2_G-pp2_G*log2(pp2_G);
    end
    if pp2_B~=0
        xxs2_B=xxs2_B-pp2_B*log2(pp2_B);
    end
end

%% 加密图像相邻图像相关性分析
%{
先随机在0~M-1行和0~N-1列选中1000个像素点，
计算水平相关性时，选择每个点的相邻的右边的点；
计算垂直相关性时，选择每个点的相邻的下方的点；
计算对角线相关性时，选择每个点的相邻的右下方的点。
%}
%相关性曲线
%水平
XX_R_SP=zeros(1,1000);YY_R_SP=zeros(1,1000);  %预分配内存
XX_G_SP=zeros(1,1000);YY_G_SP=zeros(1,1000);
XX_B_SP=zeros(1,1000);YY_B_SP=zeros(1,1000);
%垂直
XX_R_CZ=zeros(1,1000);YY_R_CZ=zeros(1,1000);  %预分配内存
XX_G_CZ=zeros(1,1000);YY_G_CZ=zeros(1,1000);
XX_B_CZ=zeros(1,1000);YY_B_CZ=zeros(1,1000);
%对角线
XX_R_DJX=zeros(1,1000);YY_R_DJX=zeros(1,1000);  %预分配内存
XX_G_DJX=zeros(1,1000);YY_G_DJX=zeros(1,1000);
XX_B_DJX=zeros(1,1000);YY_B_DJX=zeros(1,1000);
for i=1:1000
    %水平
    XX_R_SP(i)=Q_R(x1(i),y1(i));
    YY_R_SP(i)=Q_R(x1(i)+1,y1(i));
    XX_G_SP(i)=Q_G(x1(i),y1(i));
    YY_G_SP(i)=Q_G(x1(i)+1,y1(i));
    XX_B_SP(i)=Q_B(x1(i),y1(i));
    YY_B_SP(i)=Q_B(x1(i)+1,y1(i));
    %垂直
    XX_R_CZ(i)=Q_R(x1(i),y1(i));
    YY_R_CZ(i)=Q_R(x1(i),y1(i)+1);
    XX_G_CZ(i)=Q_G(x1(i),y1(i));
    YY_G_CZ(i)=Q_G(x1(i),y1(i)+1);
    XX_B_CZ(i)=Q_B(x1(i),y1(i));
    YY_B_CZ(i)=Q_B(x1(i),y1(i)+1);
    %对角线
    XX_R_DJX(i)=Q_R(x1(i),y1(i));
    YY_R_DJX(i)=Q_R(x1(i)+1,y1(i)+1);
    XX_G_DJX(i)=Q_G(x1(i),y1(i));
    YY_G_DJX(i)=Q_G(x1(i)+1,y1(i)+1);
    XX_B_DJX(i)=Q_B(x1(i),y1(i));
    YY_B_DJX(i)=Q_B(x1(i)+1,y1(i)+1);
end
%水平
figure;scatter(XX_R_SP,YY_R_SP,18,'filled');xlabel('R通道随机点像素灰度值');ylabel('与该点相邻水平方向像素灰度值');title('加密图像R通道水平相关性曲线');
figure;scatter(XX_G_SP,YY_G_SP,18,'filled');xlabel('G通道随机点像素灰度值');ylabel('与该点相邻水平方向像素灰度值');title('加密图像G通道水平相关性曲线');
figure;scatter(XX_B_SP,YY_B_SP,18,'filled');xlabel('B通道随机点像素灰度值');ylabel('与该点相邻水平方向像素灰度值');title('加密图像B通道水平相关性曲线');
%垂直
figure;scatter(XX_R_CZ,YY_R_CZ,18,'filled');xlabel('R通道随机点像素灰度值');ylabel('与该点相邻垂直方向像素灰度值');title('加密图像R通道垂直相关性曲线');
figure;scatter(XX_G_CZ,YY_G_CZ,18,'filled');xlabel('G通道随机点像素灰度值');ylabel('与该点相邻垂直方向像素灰度值');title('加密图像G通道垂直相关性曲线');
figure;scatter(XX_B_CZ,YY_B_CZ,18,'filled');xlabel('B通道随机点像素灰度值');ylabel('与该点相邻垂直方向像素灰度值');title('加密图像B通道垂直相关性曲线');
%对角线
figure;scatter(XX_R_DJX,YY_R_DJX,18,'filled');xlabel('R通道随机点像素灰度值');ylabel('与该点相邻对角线方向像素灰度值');title('加密图像R通道对角线相关性曲线');
figure;scatter(XX_G_DJX,YY_G_DJX,18,'filled');xlabel('G通道随机点像素灰度值');ylabel('与该点相邻对角线方向像素灰度值');title('加密图像G通道对角线相关性曲线');
figure;scatter(XX_B_DJX,YY_B_DJX,18,'filled');xlabel('B通道随机点像素灰度值');ylabel('与该点相邻对角线方向像素灰度值');title('加密图像B通道对角线相关性曲线');
%R通道
Q_R=double(Q_R);
EX2_R=0;EY2_SP_R=0;DX2_R=0;DY2_SP_R=0;COVXY2_SP_R=0;    %水平
EY2_CZ_R=0;DY2_CZ_R=0;COVXY2_CZ_R=0;    %垂直
EY2_DJX_R=0;DY2_DJX_R=0;COVXY2_DJX_R=0;   %对角线
%G通道
Q_G=double(Q_G);
EX2_G=0;EY2_SP_G=0;DX2_G=0;DY2_SP_G=0;COVXY2_SP_G=0;    %水平
EY2_CZ_G=0;DY2_CZ_G=0;COVXY2_CZ_G=0;    %垂直
EY2_DJX_G=0;DY2_DJX_G=0;COVXY2_DJX_G=0;   %对角线
%B通道
Q_B=double(Q_B);
EX2_B=0;EY2_SP_B=0;DX2_B=0;DY2_SP_B=0;COVXY2_SP_B=0;    %水平
EY2_CZ_B=0;DY2_CZ_B=0;COVXY2_CZ_B=0;    %垂直
EY2_DJX_B=0;DY2_DJX_B=0;COVXY2_DJX_B=0;   %对角线
for i=1:NN
    %第一个像素点的E，水平、垂直、对角线时计算得出的第一个像素点的E相同，统一用EX2表示
    EX2_R=EX2_R+Q_R(x1(i),y1(i));
    EX2_G=EX2_G+Q_G(x1(i),y1(i));
    EX2_B=EX2_B+Q_B(x1(i),y1(i));
    %第二个像素点的E，水平、垂直、对角线的E分别对应EY2_SP、EY2_CZ、EY2_DJX
    %R通道
    EY2_SP_R=EY2_SP_R+Q_R(x1(i),y1(i)+1);
    EY2_CZ_R=EY2_CZ_R+Q_R(x1(i)+1,y1(i));
    EY2_DJX_R=EY2_DJX_R+Q_R(x1(i)+1,y1(i)+1);
    %G通道
    EY2_SP_G=EY2_SP_G+Q_G(x1(i),y1(i)+1);
    EY2_CZ_G=EY2_CZ_G+Q_G(x1(i)+1,y1(i));
    EY2_DJX_G=EY2_DJX_G+Q_G(x1(i)+1,y1(i)+1);
    %B通道
    EY2_SP_B=EY2_SP_B+Q_B(x1(i),y1(i)+1);
    EY2_CZ_B=EY2_CZ_B+Q_B(x1(i)+1,y1(i));
    EY2_DJX_B=EY2_DJX_B+Q_B(x1(i)+1,y1(i)+1);
end
%统一在循环外除以像素点对数1000，可减少运算次数
%R通道
EX2_R=EX2_R/NN;
EY2_SP_R=EY2_SP_R/NN;
EY2_CZ_R=EY2_CZ_R/NN;
EY2_DJX_R=EY2_DJX_R/NN;
%G通道
EX2_G=EX2_G/NN;
EY2_SP_G=EY2_SP_G/NN;
EY2_CZ_G=EY2_CZ_G/NN;
EY2_DJX_G=EY2_DJX_G/NN;
%B通道
EX2_B=EX2_B/NN;
EY2_SP_B=EY2_SP_B/NN;
EY2_CZ_B=EY2_CZ_B/NN;
EY2_DJX_B=EY2_DJX_B/NN;

for i=1:NN
    %第一个像素点的D，水平、垂直、对角线时计算得出第一个像素点的D相同，统一用DX2表示
    DX2_R=DX2_R+(Q_R(x1(i),y1(i))-EX2_R)^2;
    DX2_G=DX2_G+(Q_G(x1(i),y1(i))-EX2_G)^2;
    DX2_B=DX2_B+(Q_B(x1(i),y1(i))-EX2_B)^2;
    %第二个像素点的D，水平、垂直、对角线的D分别对应DY2_SP、DY2_CZ、DY2_DJX
    %R通道
    DY2_SP_R=DY2_SP_R+(Q_R(x1(i),y1(i)+1)-EY2_SP_R)^2;
    DY2_CZ_R=DY2_CZ_R+(Q_R(x1(i)+1,y1(i))-EY2_CZ_R)^2;
    DY2_DJX_R=DY2_DJX_R+(Q_R(x1(i)+1,y1(i)+1)-EY2_DJX_R)^2;
    %G通道
    DY2_SP_G=DY2_SP_G+(Q_G(x1(i),y1(i)+1)-EY2_SP_G)^2;
    DY2_CZ_G=DY2_CZ_G+(Q_G(x1(i)+1,y1(i))-EY2_CZ_G)^2;
    DY2_DJX_G=DY2_DJX_G+(Q_G(x1(i)+1,y1(i)+1)-EY2_DJX_G)^2;
    %B通道
    DY2_SP_B=DY2_SP_B+(Q_B(x1(i),y1(i)+1)-EY2_SP_B)^2;
    DY2_CZ_B=DY2_CZ_B+(Q_B(x1(i)+1,y1(i))-EY2_CZ_B)^2;
    DY2_DJX_B=DY2_DJX_B+(Q_B(x1(i)+1,y1(i)+1)-EY2_DJX_B)^2;
    %两个相邻像素点相关函数的计算，水平、垂直、对角线
    %R通道
    COVXY2_SP_R=COVXY2_SP_R+(Q_R(x1(i),y1(i))-EX2_R)*(Q_R(x1(i),y1(i)+1)-EY2_SP_R);
    COVXY2_CZ_R=COVXY2_CZ_R+(Q_R(x1(i),y1(i))-EX2_R)*(Q_R(x1(i)+1,y1(i))-EY2_CZ_R);
    COVXY2_DJX_R=COVXY2_DJX_R+(Q_R(x1(i),y1(i))-EX2_R)*(Q_R(x1(i)+1,y1(i)+1)-EY2_DJX_R);
    %G通道
    COVXY2_SP_G=COVXY2_SP_G+(Q_G(x1(i),y1(i))-EX2_G)*(Q_G(x1(i),y1(i)+1)-EY2_SP_G);
    COVXY2_CZ_G=COVXY2_CZ_G+(Q_G(x1(i),y1(i))-EX2_G)*(Q_G(x1(i)+1,y1(i))-EY2_CZ_G);
    COVXY2_DJX_G=COVXY2_DJX_G+(Q_G(x1(i),y1(i))-EX2_G)*(Q_G(x1(i)+1,y1(i)+1)-EY2_DJX_G);
    %B通道
    COVXY2_SP_B=COVXY2_SP_B+(Q_B(x1(i),y1(i))-EX2_B)*(Q_B(x1(i),y1(i)+1)-EY2_SP_B);
    COVXY2_CZ_B=COVXY2_CZ_B+(Q_B(x1(i),y1(i))-EX2_B)*(Q_B(x1(i)+1,y1(i))-EY2_CZ_B);
    COVXY2_DJX_B=COVXY2_DJX_B+(Q_B(x1(i),y1(i))-EX2_B)*(Q_B(x1(i)+1,y1(i)+1)-EY2_DJX_B);
end
%统一在循环外除以像素点对数1000，可减少运算次数
%R通道
DX2_R=DX2_R/NN;
DY2_SP_R=DY2_SP_R/NN;
DY2_CZ_R=DY2_CZ_R/NN;
DY2_DJX_R=DY2_DJX_R/NN;
COVXY2_SP_R=COVXY2_SP_R/NN;
COVXY2_CZ_R=COVXY2_CZ_R/NN;
COVXY2_DJX_R=COVXY2_DJX_R/NN;
%G通道
DX2_G=DX2_G/NN;
DY2_SP_G=DY2_SP_G/NN;
DY2_CZ_G=DY2_CZ_G/NN;
DY2_DJX_G=DY2_DJX_G/NN;
COVXY2_SP_G=COVXY2_SP_G/NN;
COVXY2_CZ_G=COVXY2_CZ_G/NN;
COVXY2_DJX_G=COVXY2_DJX_G/NN;
%B通道
DX2_B=DX2_B/NN;
DY2_SP_B=DY2_SP_B/NN;
DY2_CZ_B=DY2_CZ_B/NN;
DY2_DJX_B=DY2_DJX_B/NN;
COVXY2_SP_B=COVXY2_SP_B/NN;
COVXY2_CZ_B=COVXY2_CZ_B/NN;
COVXY2_DJX_B=COVXY2_DJX_B/NN;
%水平、垂直、对角线的相关性
%R通道
RXY2_SP_R=COVXY2_SP_R/sqrt(DX2_R*DY2_SP_R);
RXY2_CZ_R=COVXY2_CZ_R/sqrt(DX2_R*DY2_CZ_R);
RXY2_DJX_R=COVXY2_DJX_R/sqrt(DX2_R*DY2_DJX_R);
%G通道
RXY2_SP_G=COVXY2_SP_G/sqrt(DX2_G*DY2_SP_G);
RXY2_CZ_G=COVXY2_CZ_G/sqrt(DX2_G*DY2_CZ_G);
RXY2_DJX_G=COVXY2_DJX_G/sqrt(DX2_G*DY2_DJX_G);
%B通道
RXY2_SP_B=COVXY2_SP_B/sqrt(DX2_B*DY2_SP_B);
RXY2_CZ_B=COVXY2_CZ_B/sqrt(DX2_B*DY2_CZ_B);
RXY2_DJX_B=COVXY2_DJX_B/sqrt(DX2_B*DY2_DJX_B);

%% 输出数据信息
disp('加密成功');
disp('密钥：');
disp(['密钥1：μ=',num2str(u),'     密钥2：x0=',num2str(x0),'    密钥3：x(0)=',num2str(X0),'    密钥4：y(0)=',num2str(Y0)]);
disp(['密钥5：z(0)=',num2str(Z0),'   密钥6：h(0)=',num2str(H0),'   密钥7：M1=',num2str(M1),'   密钥8：N1=',num2str(N1)]);
disp('信息熵：');
disp(['原始图片R通道信息熵=',num2str(xxs1_R),'  原始图片G通道信息熵=',num2str(xxs1_G),'  原始图片B通道信息熵=',num2str(xxs1_B)]);
disp(['加密图片R通道信息熵=',num2str(xxs2_R),'  加密图片G通道信息熵=',num2str(xxs2_G),'  加密图片B通道信息熵=',num2str(xxs2_B)]);
disp('R通道相关性：');
disp(['原始图片R通道相关性：','  水平相关性=',num2str(RXY1_SP_R),'    垂直相关性=',num2str(RXY1_CZ_R),'  对角线相关性=',num2str(RXY1_DJX_R)]);
disp(['加密图片R通道相关性：','  水平相关性=',num2str(RXY2_SP_R),'  垂直相关性=',num2str(RXY2_CZ_R),'  对角线相关性=',num2str(RXY2_DJX_R)]);
disp('G通道相关性：');
disp(['原始图片G通道相关性：','  水平相关性=',num2str(RXY1_SP_G),'   垂直相关性=',num2str(RXY1_CZ_G),'  对角线相关性=',num2str(RXY1_DJX_G)]);
disp(['加密图片G通道相关性：','  水平相关性=',num2str(RXY2_SP_G),'  垂直相关性=',num2str(RXY2_CZ_G),'  对角线相关性=',num2str(RXY2_DJX_G)]);
disp('B通道相关性：');
disp(['原始图片B通道相关性：','  水平相关性=',num2str(RXY1_SP_B),'   垂直相关性=',num2str(RXY1_CZ_B),'  对角线相关性=',num2str(RXY1_DJX_B)]);
disp(['加密图片B通道相关性：','  水平相关性=',num2str(RXY2_SP_B),'  垂直相关性=',num2str(RXY2_CZ_B),'  对角线相关性=',num2str(RXY2_DJX_B)]);