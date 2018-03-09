%% 子函数 分块函数，t为每块的边长，I为要分块的图像，num为返回第几大块
function fv=fenkuai(t,I,num)
[~,N]=size(I);
N=N/t;
x=floor(num/N)+1;      %第几大行
y=mod(num,N);           %第几大列
if y==0
    x=x-1;
    y=N;
end
fv=I(t*(x-1)+1:t*x,t*(y-1)+1:t*y);


