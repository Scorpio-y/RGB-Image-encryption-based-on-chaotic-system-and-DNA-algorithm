%% 基于混沌系统与DNA编码的彩色数字图像加密解密
clear;clc;
I=imread('../原始、加密、解密图片/加密后的lena.png','png');           %读取图像信息
I1=I(:,:,1);     %R通道
I2=I(:,:,2);     %G通道
I3=I(:,:,3);     %B通道
[M,N]=size(I1);                      %将图像的行列赋值给M,N
t=4;    %分块大小
M1=0;   %加密时补零的参数，M1=mod(M,t);作为密钥
N1=0;   %加密时补零的参数，N1=mod(N,t);作为密钥
SUM=M*N;
u=3.99;
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

for i=N:-1:1
    temp = I1(:,i);
    I1(:,i) = I1(:,Uy(i));
    I1(:,Uy(i)) = temp;
    temp = I2(:,i);
    I2(:,i) = I2(:,Uy(i));
    I2(:,Uy(i)) = temp;
    temp = I3(:,i);
    I3(:,i) = I3(:,Uy(i));
    I3(:,Uy(i)) = temp;
end
for i=M:-1:1
    temp = I1(i,:);
    I1(i,:) = I1(Ux(i),:);
    I1(Ux(i),:) = temp;
    temp = I2(i,:);
    I2(i,:) = I2(Ux(i),:);
    I2(Ux(i),:) = temp;
    temp = I3(i,:);
    I3(i,:) = I3(Ux(i),:);
    I3(Ux(i),:) = temp;
end
%% 2.产生Logistic混沌序列
% u=3.990000000000001; %密钥敏感性测试  10^-15
u=3.99;%密钥：Logistic参数μ
% x0=0.7067000000000001; %密钥敏感性测试  10^-16
x0=0.7067; %密钥：Logistic初值x0
% x0=0.3462;            %home图片
p=zeros(1,SUM+1000);
p(1)=x0;
for i=1:SUM+999                        %进行N-1次循环
    p(i+1)=u*p(i)*(1-p(i));          %循环产生密码
end
p=p(1001:length(p));

%% 3.将p序列变换到0~255范围内整数，转换成M*N的二维矩阵R
p=mod(ceil(p*10^3),256);
R=reshape(p,N,M)';  %转成M行N列

%% 4.求解混沌方程
%求四个初值X0,Y0,Z0,H0
r=(M/t)*(N/t);
% X0=0.5008000000000001;        %密钥敏感性测试
X0=0.5008;
Y0=0.5109;
Z0=0.4893;
H0=0.7765;
% X0=0.5056;        %home图片
% Y0=0.505;
% Z0=0.4564;
% H0=0.3062;
A=chen_output(X0,Y0,Z0,H0,r);
X=A(:,1);
X=X(1502:length(X));
Y=A(:,2);
Y=Y(1502:length(Y));
Z=A(:,3);
Z=Z(1502:length(Z));
H=A(:,4);
H=H(1502:length(H));

%% 5.DNA编码
%X,Y分别决定I和R的DNA编码方式，有8种，1~8
X=mod(floor(X*10^4),8)+1;
Y=mod(floor(Y*10^4),8)+1;
Z=mod(floor(Z*10^4),3);
Z(Z==0)=3;
Z(Z==1)=0;
Z(Z==3)=1;
H=mod(floor(H*10^4),8)+1;
e=N/t;
for i=r:-1:2
    Q1_R=DNA_bian(fenkuai(t,I1,i),H(i));
    Q1_G=DNA_bian(fenkuai(t,I2,i),H(i));
    Q1_B=DNA_bian(fenkuai(t,I3,i),H(i));
    
    Q1_last_R=DNA_bian(fenkuai(t,I1,i-1),H(i-1));
    Q1_last_G=DNA_bian(fenkuai(t,I2,i-1),H(i-1));
    Q1_last_B=DNA_bian(fenkuai(t,I3,i-1),H(i-1));
    
    Q2_R=DNA_yunsuan(Q1_R,Q1_last_R,Z(i));        %扩散前
    Q2_G=DNA_yunsuan(Q1_G,Q1_last_G,Z(i));
    Q2_B=DNA_yunsuan(Q1_B,Q1_last_B,Z(i));

    Q3=DNA_bian(fenkuai(t,R,i),Y(i));
    
    Q4_R=DNA_yunsuan(Q2_R,Q3,Z(i));
    Q4_G=DNA_yunsuan(Q2_G,Q3,Z(i));
    Q4_B=DNA_yunsuan(Q2_B,Q3,Z(i));
    
    xx=floor(i/e)+1;
    yy=mod(i,e);
    if yy==0
        xx=xx-1;
        yy=e;
    end
    I1((xx-1)*t+1:xx*t,(yy-1)*t+1:yy*t)=DNA_jie(Q4_R,X(i));
    I2((xx-1)*t+1:xx*t,(yy-1)*t+1:yy*t)=DNA_jie(Q4_G,X(i));
    I3((xx-1)*t+1:xx*t,(yy-1)*t+1:yy*t)=DNA_jie(Q4_B,X(i));
end
Q5_R=DNA_bian(fenkuai(t,I1,1),H(1));
Q5_G=DNA_bian(fenkuai(t,I2,1),H(1));
Q5_B=DNA_bian(fenkuai(t,I3,1),H(1));

Q6=DNA_bian(fenkuai(t,R,1),Y(1));

Q7_R=DNA_yunsuan(Q5_R,Q6,Z(1));
Q7_G=DNA_yunsuan(Q5_G,Q6,Z(1));
Q7_B=DNA_yunsuan(Q5_B,Q6,Z(1));

I1(1:t,1:t)=DNA_jie(Q7_R,X(1));
I2(1:t,1:t)=DNA_jie(Q7_G,X(1));
I3(1:t,1:t)=DNA_jie(Q7_B,X(1));

Q_jiemi(:,:,1)=uint8(I1);
Q_jiemi(:,:,2)=uint8(I2);
Q_jiemi(:,:,3)=uint8(I3);

%% 6、去除加密时补的零
if M1~=0
    Q_jiemi=Q_jiemi(1:M-t+M1,:,:);
end
if N1~=0
    Q_jiemi=Q_jiemi(:,1:N-t+N1,:);
end
%比较解密后的图与原图是否完全相同
II=imread('../原始、加密、解密图片/lena.png','png');
cha=sum(sum(sum(Q_jiemi-II)));      %两幅图做差后求总和
%% 保存图片
imwrite(Q_jiemi,'../原始、加密、解密图片/解密后的lena.png','png');       
disp('您输入的解密密钥为：');
disp(['密钥1：μ=',num2str(u),'     密钥2：x0=',num2str(x0),'    密钥3：x(0)=',num2str(X0),'    密钥4：y(0)=',num2str(Y0)]);
disp(['密钥5：z(0)=',num2str(Z0),'   密钥6：h(0)=',num2str(H0),'   密钥7：M1=',num2str(M1),'   密钥8：N1=',num2str(N1)]);
disp('解密完成'); 
imshow(Q_jiemi);title('解密后图片');