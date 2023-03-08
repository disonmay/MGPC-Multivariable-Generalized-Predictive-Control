function y=mconv(x1,x2)
m=length(x1(1,:))/2;
n=length(x2(1,:))/2;
y=zeros(2,2);
for j=1:n
    for i=1:m
        x3(:,2*(i-1)+1:2*i)=x1(:,2*(i-1)+1:2*i)*x2(:,2*(j-1)+1:2*j);           
    end   
    temp=[zeros(2,2*(j-1)),x3]; 
    y=polyadd(y,temp);
end