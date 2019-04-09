function f=chao_SimpleLorenz(t,y)

f=zeros(4,1);
f(1)=35*(y(2)-y(1))+y(4);
f(2)=7*y(1)-y(1)*y(3)+12*y(2);
f(3)=y(1)*y(2)-3*y(3);
f(4)=y(2)*y(3)+0.58*y(4);