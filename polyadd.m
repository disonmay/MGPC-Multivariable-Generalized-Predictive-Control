function y=polyadd(x1,x2)
n1=length(x1(1,:));
n2=length(x2(1,:));
if n1>n2
    x2=[x2,zeros(2,n1-n2)];
elseif n1<n2
    x1=[x1,zeros(2,n2-n1)];
end
y=x1+x2;
