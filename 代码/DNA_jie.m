%% 子函数 DNA解码
function fv=DNA_jie(array,num)
[m,n]=size(array);
A=zeros(m,n);           %预分配内存
if num==1 || num==2
    for i=1:m
        for j=1:n
            if array(i,j)=='A'
                A(i,j)=0;
            elseif array(i,j)=='T'
                A(i,j)=3;
            elseif array(i,j)=='G'
                if num==1
                    A(i,j)=1;
                else
                    A(i,j)=2;
                end
            else
                if num==1
                    A(i,j)=2;
                else
                    A(i,j)=1;
                end
            end
        end
    end
elseif num==3 || num==4
    for i=1:m
        for j=1:n
            if array(i,j)=='A'
                A(i,j)=1;
            elseif array(i,j)=='T'
                A(i,j)=2;
            elseif array(i,j)=='G'
                if num==3
                    A(i,j)=0;
                else
                    A(i,j)=3;
                end
            else
                if num==3
                    A(i,j)=3;
                else
                    A(i,j)=0;
                end
            end
        end
    end
elseif num==5 || num==6
    for i=1:m
        for j=1:n
            if array(i,j)=='A'
                A(i,j)=2;
            elseif array(i,j)=='T'
                A(i,j)=1;
            elseif array(i,j)=='G'
                if num==5
                    A(i,j)=0;
                else
                    A(i,j)=3;
                end
            else
                if num==5
                    A(i,j)=3;
                else
                    A(i,j)=0;
                end
            end
        end
    end
else                    %  num==7 || num==8
    for i=1:m
        for j=1:n
            if array(i,j)=='A'
                A(i,j)=3;
            elseif array(i,j)=='T'
                A(i,j)=0;
            elseif array(i,j)=='G'
                if num==7
                    A(i,j)=1;
                else
                    A(i,j)=2;
                end
            else
                if num==7
                    A(i,j)=2;
                else
                    A(i,j)=1;
                end
            end
        end
    end
end
A1=A(:,1:m);
A2=A(:,m+1:2*m);
A3=A(:,2*m+1:3*m);
A4=A(:,3*m+1:4*m);
fv=A1*64+A2*16+A3*4+A4;
end