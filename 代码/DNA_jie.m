%% ×Óº¯Êı DNA½âÂë
function fv=DNA_jie(array,num)
[m,n]=size(array);
if num==1
    for i=1:m
        for j=1:n
            if array(i,j)=='A'
                A(i,j)=0;
            elseif array(i,j)=='T'
                A(i,j)=3;
            elseif array(i,j)=='G'
                A(i,j)=2;
            else
                A(i,j)=1;
            end
        end
    end
elseif num==2
    for i=1:m
        for j=1:n
            if array(i,j)=='A'
                A(i,j)=0;
            elseif array(i,j)=='T'
                A(i,j)=3;
            elseif array(i,j)=='G'
                A(i,j)=1;
            else
                A(i,j)=2;
            end
        end
    end
elseif num==3
    for i=1:m
        for j=1:n
            if array(i,j)=='A'
                A(i,j)=3;
            elseif array(i,j)=='T'
                A(i,j)=0;
            elseif array(i,j)=='G'
                A(i,j)=2;
            else
                A(i,j)=1;
            end
        end
    end
elseif num==4
    for i=1:m
        for j=1:n
            if array(i,j)=='A'
                A(i,j)=3;
            elseif array(i,j)=='T'
                A(i,j)=0;
            elseif array(i,j)=='G'
                A(i,j)=1;
            else
                A(i,j)=2;
            end
        end
    end
elseif num==5
    for i=1:m
        for j=1:n
            if array(i,j)=='A'
                A(i,j)=1;
            elseif array(i,j)=='T'
                A(i,j)=2;
            elseif array(i,j)=='G'
                A(i,j)=0;
            else
                A(i,j)=3;
            end
        end
    end
elseif num==6
    for i=1:m
        for j=1:n
            if array(i,j)=='A'
                A(i,j)=2;
            elseif array(i,j)=='T'
                A(i,j)=1;
            elseif array(i,j)=='G'
                A(i,j)=0;
            else
                A(i,j)=3;
            end
        end
    end
elseif num==7
    for i=1:m
        for j=1:n
            if array(i,j)=='A'
                A(i,j)=1;
            elseif array(i,j)=='T'
                A(i,j)=2;
            elseif array(i,j)=='G'
                A(i,j)=3;
            else
                A(i,j)=0;
            end
        end
    end
else
    for i=1:m
        for j=1:n
            if array(i,j)=='A'
                A(i,j)=2;
            elseif array(i,j)=='T'
                A(i,j)=1;
            elseif array(i,j)=='G'
                A(i,j)=3;
            else
                A(i,j)=0;
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