%% 子函数 DNA运算 arr1和arr2表示两个要运算的DNA矩阵。
%num为0时，表示加；num为1时，表示减；num为2时，表示异或
function fv=DNA_yunsuan(arr1,arr2,num)
[m,n]=size(arr1);
if num==0    %加法
    for i=1:m
        for j=1:n
            if arr1(i,j)=='A'
                fv(i,j)=arr2(i,j);
            elseif arr1(i,j)=='T'
                if arr2(i,j)=='A'
                    fv(i,j)='T';
                elseif arr2(i,j)=='T'
                    fv(i,j)='C';
                elseif arr2(i,j)=='G'
                    fv(i,j)='A';
                else
                    fv(i,j)='G';
                end
            elseif arr1(i,j)=='G'
                if arr2(i,j)=='A'
                    fv(i,j)='G';
                elseif arr2(i,j)=='T'
                    fv(i,j)='A';
                elseif arr2(i,j)=='G'
                    fv(i,j)='C';
                else
                    fv(i,j)='T';
                end
            else
                if arr2(i,j)=='A'
                    fv(i,j)='C';
                elseif arr2(i,j)=='T'
                    fv(i,j)='G';
                elseif arr2(i,j)=='G'
                    fv(i,j)='T';
                else
                    fv(i,j)='A';
                end
            end
        end
    end
elseif num==1       %减法
    for i=1:m
        for j=1:n
            if arr1(i,j)=='A'
                if arr2(i,j)=='A'
                    fv(i,j)='A';
                elseif arr2(i,j)=='T'
                    fv(i,j)='G';
                elseif arr2(i,j)=='G'
                    fv(i,j)='T';
                else
                    fv(i,j)='C';
                end
            elseif arr1(i,j)=='T'
                if arr2(i,j)=='A'
                    fv(i,j)='T';
                elseif arr2(i,j)=='T'
                    fv(i,j)='A';
                elseif arr2(i,j)=='G'
                    fv(i,j)='C';
                else
                    fv(i,j)='G';
                end
            elseif arr1(i,j)=='G'
                if arr2(i,j)=='A'
                    fv(i,j)='G';
                elseif arr2(i,j)=='T'
                    fv(i,j)='C';
                elseif arr2(i,j)=='G'
                    fv(i,j)='A';
                else
                    fv(i,j)='T';
                end
            else
                if arr2(i,j)=='A'
                    fv(i,j)='C';
                elseif arr2(i,j)=='T'
                    fv(i,j)='T';
                elseif arr2(i,j)=='G'
                    fv(i,j)='G';
                else
                    fv(i,j)='A';
                end
            end
        end
    end
else        %异或
    for i=1:m
        for j=1:n
            if arr1(i,j)=='A'
                    fv(i,j)=arr2(i,j);
            elseif arr1(i,j)=='T'
                if arr2(i,j)=='A'
                    fv(i,j)='T';
                elseif arr2(i,j)=='T'
                    fv(i,j)='A';
                elseif arr2(i,j)=='G'
                    fv(i,j)='C';
                else
                    fv(i,j)='G';
                end
            elseif arr1(i,j)=='G'
                if arr2(i,j)=='A'
                    fv(i,j)='G';
                elseif arr2(i,j)=='T'
                    fv(i,j)='C';
                elseif arr2(i,j)=='G'
                    fv(i,j)='A';
                else
                    fv(i,j)='T';
                end
            else
                if arr2(i,j)=='A'
                    fv(i,j)='C';
                elseif arr2(i,j)=='T'
                    fv(i,j)='G';
                elseif arr2(i,j)=='G'
                    fv(i,j)='T';
                else
                    fv(i,j)='A';
                end
            end
        end
    end
end
