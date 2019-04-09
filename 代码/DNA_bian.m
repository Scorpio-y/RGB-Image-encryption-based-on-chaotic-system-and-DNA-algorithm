%% 子函数 DNA编码
function fv=DNA_bian(array,num)
t=length(array);
a1=bitand(array,192)/64;
a2=bitand(array,48)/16;
a3=bitand(array,12)/4;
a4=bitand(array,3);
A=[a1,a2,a3,a4];        % t行4t列
fv=zeros(t,4*t);        % 预分配内存
if num==1 || num==2
    for i=1:t
        for j=1:4*t
            if A(i,j)==0
                fv(i,j)='A';
            elseif A(i,j)==3
                fv(i,j)='T';
            elseif A(i,j)==2
                if num==1
                    fv(i,j)='C';
                else
                    fv(i,j)='G';
                end
            else
                if num==1
                    fv(i,j)='G';
                else
                    fv(i,j)='C';
                end
            end
        end
    end
elseif num==3 || num==4
    for i=1:t
        for j=1:4*t
            if A(i,j)==1
                fv(i,j)='A';
            elseif A(i,j)==2
                fv(i,j)='T';
            elseif A(i,j)==0
                if num==3
                    fv(i,j)='G';
                else
                    fv(i,j)='C';
                end
            else
                if num==3
                    fv(i,j)='C';
                else
                    fv(i,j)='G';
                end
            end
        end
    end
elseif num==5 || num==6
    for i=1:t
        for j=1:4*t
            if A(i,j)==2
                fv(i,j)='A';
            elseif A(i,j)==1
                fv(i,j)='T';
            elseif A(i,j)==0
                if num==5
                    fv(i,j)='G';
                else
                    fv(i,j)='C';
                end
            else
                if num==5
                    fv(i,j)='C';
                else
                    fv(i,j)='G';
                end
            end
        end
    end
else                    %  num==7 || num==8
    for i=1:t
        for j=1:4*t
            if A(i,j)==3
                fv(i,j)='A';
            elseif A(i,j)==0
                fv(i,j)='T';
            elseif A(i,j)==1
                if num==7
                    fv(i,j)='G';
                else
                    fv(i,j)='C';
                end
            else
                if num==7
                    fv(i,j)='C';
                else
                    fv(i,j)='G';
                end
            end
        end
    end
end