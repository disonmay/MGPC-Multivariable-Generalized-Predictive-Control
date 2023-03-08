%% 自校正带约束MGPC
clear;clc;
%% 数据初始化，代入题干
yr=[1;1];
arfa=0.2;
k=1;
epoch=100;
na=2;nb=2;
A0=[1,0;0,1];
A1=[-1.4,-0.2;-0.1,-0.9];
A2=[0.48,0.1;0,0.2];
B0=[1,0;0,1];
B1=[1.5,1;0,1];
AAA=[A0,A1,A2];
BBB=[B0,B1];
A0=[1,0;0,1];
A1=[ones(2,2)/1000];
A2=[ones(2,2)/1000];
A=[A0,A1,A2];
B0=[ones(2,2)/1000];
B1=[ones(2,2)/1000];
B=[B0,B1];
N=input('请输入N:');
Nu=input('请输入Nu:');
lamda=input('请输入lamda:');
detau_ys=input('请输入Δu(k)的约束条件:'); % 建议输入[1 1]
u_ys=input('请输入u(k)的约束条件:'); % 建议输入[1.5 1.5]
y=[zeros(1,2),zeros(1,2*na)]';u=zeros(1,2*nb)';
c0=[A(:,3:2*(na+1)),B]';
p0=10^6*eye(2*na+2*nb,2*na+2*nb);
yy=[zeros(1,2),zeros(1,2*na)]';
drtuu=zeros(1,2*nb)';
%% 求解输出信号
for i=1:k:epoch
    kesi=-0.05+0.1*rand(2,1);
    yk=(-AAA)*y+BBB*u+kesi;
    yna=yy(2*na+1:2*na+2)';
    yy=[yk',y(3:2*(na+1))']';
    yyy=[y(3:2*3)',yna]';
    drtyy=yy-yyy;
    for l=na:-1:2
        y(2*l+1:2*l+2,1)=y(2*(l-1)+1:2*(l-1)+2,1);
    end
    y(3:4,1)=yk;
    yd=yk;
    yK(:,i)=yk;
    % 参数辨识
    for j=1:2
        h1=[-drtyy(3:2*3)',drtuu']'; %辨识
        x=h1'*p0*h1+1; 
        x1=inv(x);
        k1=p0*h1*x1;
        d1(j)=drtyy(j)-h1'*c0(:,j); 
        c1(:,j)=c0(:,j)+k1*d1(j);
        c0(:,j)=c1(:,j);
        p1=p0-k1*k1'*(h1'*p0*h1+1);
        p0=p1;
    end
        A=[eye(2,2),c1(1:2*na,:)'];B=c1(2*na+1:2*na+2*nb,:)';
%% 丢番图方程组的递推
ta=[1 0 -1 0
    0 1 0 -1];%求G
AA=mconv(A,ta);
E=[1 0
    0 1];
EE=[E,zeros(2,2*(N-1))];
F=-AA(:,3:end);
FF=F;
EB=mconv(E,B);
G=[fliplr(EB(:,1:2)),zeros(2,2*(Nu-1))];
H=EB(:,3:end);
for j=2:1:N
    zj=[zeros(2,2*(j-1)),eye(2,2)];
    ej=F(:,1:2);
    E=polyadd(E,mconv(ej,zj));
    f=polyadd(F,mconv(-ej,AA));
    F=f(:,3:end);
    EE(2*(j-1)+1:2*(j-1)+2,:)=[E,zeros(2,2*(N-j))];
    FF(2*(j-1)+1:2*(j-1)+2,:)=F;
    EB=mconv(E,B);
    if j<=Nu
        G(2*(j-1)+1:2*(j-1)+2,:)=[fliplr(EB(:,1:2*j)),zeros(2,2*(Nu-j))];  
    elseif j>Nu
        G(2*(j-1)+1:2*(j-1)+2,:)=fliplr(EB(:,2*(j-Nu+1)-1:2*j));
    end
    H(2*(j-1)+1:2*(j-1)+2,:)=EB(:,2*j+1:end);
end
for l=1:2:2*Nu-1
    G(:,l:l+1)=fliplr(G(:,l:l+1));
end
GG=inv(G'*G+lamda*eye(2*Nu,2*Nu))*G';
    % 柔化因子α柔化序列向量
    for j=1:1:N
        yd=arfa*yd+(1-arfa)*yr;
        ydk(2*j-1:2*j,1)=yd;
        y0=FF(2*(j-1)+1:2*(j-1)+2,:)*yy+H(2*(j-1)+1:2*(j-1)+2,:)*(u(1:2,1)-u(3:4,1));
        y0k(2*j-1:2*j,1)=y0;
    end
    drtu=GG*(ydk-y0k);
    for e=1:2
        if drtu(e,1)>detau_ys(1,e)
            drtu(e,1)=detau_ys(1,e);
        elseif drtu(e,1)<(-1*detau_ys(1,e))
            drtu(e,1)=-1*detau_ys(1,e);
        end
    end
    u0=u;
    for m=nb-1:-1:1
        u(2*m+1:2*m+2,1)=u(2*(m-1)+1:2*(m-1)+2,1);
    end
    u(1:2,1)=drtu(1:2,1)+u(1:2,1);
    for e=1:2
        if u(e,1)>u_ys(1,e)
            u(e,1)=u_ys(1,e);
        elseif u(e,1)<(-1*u_ys(1,e))
            u(e,1)=-1*u_ys(1,e);
        end
    end
    drtuu=u-u0;
    U(:,i)=u(1:2,1);
end
%% 曲线展示
i=1:k:epoch; 
yR=zeros(2,length(i));
yR(1,:)=ones(1,length(i)).*yr(1);
yR(2,:)=ones(1,length(i)).*yr(2);
subplot(2,2,1);
 plot(i,yR(1,:),'r',i,yK(1,:),'g');
 xlabel('k');
ylabel('y1(k)');
legend('设定值','实际追踪值');
title('y1');
subplot(2,2,2);
 plot(i,yR(2,:),'r',i,yK(2,:),'g');
 xlabel('k');
ylabel('y2(k)');
legend('设定值','实际追踪值');
title('y2');
subplot(2,2,3);
 plot(i,U(1,:));
 xlabel('k');
ylabel('u1(k)');
title('u1');
subplot(2,2,4);
 plot(i,U(2,:));
 xlabel('k');
ylabel('u2(k)');
title('u2');