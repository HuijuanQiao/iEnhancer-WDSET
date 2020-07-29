function [x,k]=DS_fusion(x,y)
[nx,mx]=size(x);
if 1~=nx
    disp('x must be row vector');
    return;
end
[ny,my]=size(y);
if 1~=ny
     disp('y must be row vector');
    return;
end
if mx~=my
    disp('The column of x and y must be equal');
    return;
end
temp=0;
for i=1:mx-1
    
    if i==mx-1
        x(1,i)=x(1,i)*y(1,i); 
    else
        x(1,i)=x(1,i)*y(1,i)+x(1,i)*y(1,mx-1)+y(1,i)*x(1,mx-1);
    end
    temp=temp+x(1,i);
end
for i=1:mx-1
    x(1,i)=x(1,i)/temp;
end
x(1,mx)=0;
k=1-temp;
end