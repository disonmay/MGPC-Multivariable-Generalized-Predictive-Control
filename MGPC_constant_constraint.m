%% 定常带约束MGPC
clear;clc;
%% 数据初始化，代入题干
arfa=0.2; % 柔化因子α
k=1; % 仿真间隔
epoch=100; % 仿真步长
w=input("请输入给定信号类型：1——方波；2——阶跃；3——正弦\n");
switch w
    case 1
        yr(1,:)=5*[ones(epoch/4,1);-ones(epoch/4,1);ones(epoch/4,1);-ones(epoch/4+1,1)];
        yr(2,:)=10*[ones(epoch/4,1);-ones(epoch/4,1);ones(epoch/4,1);-ones(epoch/4+1,1)];
    case 2
        yr(1,:)=5*[ones(epoch/4,1);ones(epoch/4,1);ones(epoch/4,1);ones(epoch/4+1,1)];
        yr(2,:)=10*[ones(epoch/4,1);ones(epoch/4,1);ones(epoch/4,1);ones(epoch/4+1,1)];
    case 3
        for i=1:epoch+1                                           
            p(i)=5*sin(2*pi*i/20);
            yr=p;
        end
end
na=2;nb=2; % A和B的阶次
A0=[1,0;0,1]; % A0的系数矩阵
A1=[-1.4,-0.2;-0.1,-0.9]; % A1的系数矩阵
A2=[0.48,0.1;0,0.2]; %A2的系数矩阵
B0=[1,0;0,1]; % B0的系数矩阵
B1=[1.5,1;0,1]; % B1的系数矩阵
A=[A0,A1,A2]; % A的系数矩阵化为多项式
B=[B0,B1]; % B的系数矩阵化为多项式
% 输入三个初始数据：预测长度、控制长度和lamda
N=input('请输入N:');
Nu=input('请输入Nu:');
lamda=input('请输入lamda:');
detau_ys=input('请输入Δu(k)的约束条件:'); % 建议输入[1 1]
u_ys=input('请输入u(k)的约束条件:'); % 建议输入[3 3]
y=[zeros(1,2),zeros(1,2*na)]'; % y=[0 0 0 0 0 0]'
u=zeros(1,2*nb)'; % u=[0 0 0 0]'——控制信号，前两项u(1:2,1)为u(k-1)，后两项u(3:4,1)为u(k-2)
%% 丢番图的递推
ta=[1 0 -1 0
    0 1 0 -1];
AA=mconv(A,ta); % mconv为多项式乘积函数
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
for i=1:2:2*Nu-1
    G(:,i:i+1)=fliplr(G(:,i:i+1));
end
GG=inv(G'*G+lamda*eye(2*Nu,2*Nu))*G'; % 求解ΔU(k)时的系数多项式
%% 求解输出与控制信号
for i=1:k:epoch
    kesi=-0.05+0.1*rand(2,1); % 噪声向量ξ
    yk=(-A)*y+B*u+kesi; % 被控对象模型，同时求解出k时刻即初始时的实际输出
    yy=[yk',y(3:2*(na+1))']'; % yy=[y(k),y(k-1),y(k-2)]
    for l=na:-1:2
        y(2*l+1:2*l+2,1)=y(2*(l-1)+1:2*(l-1)+2,1); % y(5:6,1)=y(3:4,1)——>y(k-1)的值赋给y(k-2)
    end
    y(3:4,1)=yk; % 此时刻y(k)成为下一时刻y(k-1)
    % 柔化因子α柔化序列向量
    yd=yk;
    yK(:,i)=yk; % 存储每步的输出
    for j=1:1:N
        yd=arfa*yd+(1-arfa)*yr(:,i+1);
        ydk(2*j-1:2*j,1)=yd;
        y0=FF(2*(j-1)+1:2*(j-1)+2,:)*yy+H(2*(j-1)+1:2*(j-1)+2,:)*(u(1:2,1)-u(3:4,1));
        y0k(2*j-1:2*j,1)=y0;
    end
    % 计算ΔU(k)
    drtu=GG*(ydk-y0k);
    for e=1:2
        if drtu(e,1)>detau_ys(1,e)
            drtu(e,1)=detau_ys(1,e);
        elseif drtu(e,1)<(-1*detau_ys(1,e))
            drtu(e,1)=-1*detau_ys(1,e);
        end
    end
    for m=nb-1:-1:1
        u(2*m+1:2*m+2,1)=u(2*(m-1)+1:2*(m-1)+2,1); % u(3:4,1)=u(1:2,1)
    end
    u(1:2,1)=drtu(1:2,1)+u(1:2,1); % u(1:2,1)=Δu(1:2,1)+u(1:2,1)
    for e=1:2
        if u(e,1)>u_ys(1,e)
            u(e,1)=u_ys(1,e);
        elseif u(e,1)<(-1*u_ys(1,e))
            u(e,1)=-1*u_ys(1,e);
        end
    end
    U(:,i)=u(1:2,1); % 储存控制信号
end
%% 曲线展示
i=1:k:epoch; %画图
if w==3
    yR(1,:)=yr(1:epoch);
    yR(2,:)=yr(1:epoch);
else
    yR(1,:)=yr(1,1:epoch);
    yR(2,:)=yr(2,1:epoch);
end
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
title('u2');
ylabel('u2(k)');