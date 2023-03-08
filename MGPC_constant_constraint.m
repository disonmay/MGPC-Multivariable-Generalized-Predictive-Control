%% ������Լ��MGPC
clear;clc;
%% ���ݳ�ʼ�����������
arfa=0.2; % �ữ���Ӧ�
k=1; % ������
epoch=100; % ���沽��
w=input("����������ź����ͣ�1����������2������Ծ��3��������\n");
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
na=2;nb=2; % A��B�Ľ״�
A0=[1,0;0,1]; % A0��ϵ������
A1=[-1.4,-0.2;-0.1,-0.9]; % A1��ϵ������
A2=[0.48,0.1;0,0.2]; %A2��ϵ������
B0=[1,0;0,1]; % B0��ϵ������
B1=[1.5,1;0,1]; % B1��ϵ������
A=[A0,A1,A2]; % A��ϵ������Ϊ����ʽ
B=[B0,B1]; % B��ϵ������Ϊ����ʽ
% ����������ʼ���ݣ�Ԥ�ⳤ�ȡ����Ƴ��Ⱥ�lamda
N=input('������N:');
Nu=input('������Nu:');
lamda=input('������lamda:');
detau_ys=input('�����릤u(k)��Լ������:'); % ��������[1 1]
u_ys=input('������u(k)��Լ������:'); % ��������[3 3]
y=[zeros(1,2),zeros(1,2*na)]'; % y=[0 0 0 0 0 0]'
u=zeros(1,2*nb)'; % u=[0 0 0 0]'���������źţ�ǰ����u(1:2,1)Ϊu(k-1)��������u(3:4,1)Ϊu(k-2)
%% ����ͼ�ĵ���
ta=[1 0 -1 0
    0 1 0 -1];
AA=mconv(A,ta); % mconvΪ����ʽ�˻�����
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
GG=inv(G'*G+lamda*eye(2*Nu,2*Nu))*G'; % ��⦤U(k)ʱ��ϵ������ʽ
%% ������������ź�
for i=1:k:epoch
    kesi=-0.05+0.1*rand(2,1); % ����������
    yk=(-A)*y+B*u+kesi; % ���ض���ģ�ͣ�ͬʱ����kʱ�̼���ʼʱ��ʵ�����
    yy=[yk',y(3:2*(na+1))']'; % yy=[y(k),y(k-1),y(k-2)]
    for l=na:-1:2
        y(2*l+1:2*l+2,1)=y(2*(l-1)+1:2*(l-1)+2,1); % y(5:6,1)=y(3:4,1)����>y(k-1)��ֵ����y(k-2)
    end
    y(3:4,1)=yk; % ��ʱ��y(k)��Ϊ��һʱ��y(k-1)
    % �ữ���Ӧ��ữ��������
    yd=yk;
    yK(:,i)=yk; % �洢ÿ�������
    for j=1:1:N
        yd=arfa*yd+(1-arfa)*yr(:,i+1);
        ydk(2*j-1:2*j,1)=yd;
        y0=FF(2*(j-1)+1:2*(j-1)+2,:)*yy+H(2*(j-1)+1:2*(j-1)+2,:)*(u(1:2,1)-u(3:4,1));
        y0k(2*j-1:2*j,1)=y0;
    end
    % ���㦤U(k)
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
    u(1:2,1)=drtu(1:2,1)+u(1:2,1); % u(1:2,1)=��u(1:2,1)+u(1:2,1)
    for e=1:2
        if u(e,1)>u_ys(1,e)
            u(e,1)=u_ys(1,e);
        elseif u(e,1)<(-1*u_ys(1,e))
            u(e,1)=-1*u_ys(1,e);
        end
    end
    U(:,i)=u(1:2,1); % ��������ź�
end
%% ����չʾ
i=1:k:epoch; %��ͼ
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
legend('�趨ֵ','ʵ��׷��ֵ');
title('y1');
subplot(2,2,2);
plot(i,yR(2,:),'r',i,yK(2,:),'g');
xlabel('k');
ylabel('y2(k)');
legend('�趨ֵ','ʵ��׷��ֵ');
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