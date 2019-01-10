function [track,direction]=trajectory(v_ini,tf) 

dt=0.01;
t=0:dt:tf;
N=size(t,2);
X=zeros(1,N);
Y=X;

for i=1:N
    Z1=(2.4/25)*(X(1,i)-27.19)-1.2;
    Z2=(2.4/21.95)*(X(1,i)-56.46)-1.2;
    
    Y(1,i)=4.05/2*(1+tanh(Z1))-5.7/2*(1+tanh(Z2));
    %Y(1,i)=4.05/2*(1+tanh(2.4/25*(X(1,i)-27.19)-1.2))-5.7/2*(1+tanh(2.4/21.95*(X(1,i)-56.46)-1.2));
    dydx(1,i)=(4.05/2)*((1-tanh(Z1)^2)*2.4/25)-5.7/2*((1-tanh(Z2)^2)*2.4/21.95);
    %dydx(1,i)= (4.05/2)*(1-2.4/25*(tanh(2.4/25*(X(1,i)-27.19)-1.2))^2)-(5.7/2)*(1-2.4/21.95*(tanh(2.4/21.95*(X(1,i)-56.46)-1.2))^2);
    
    dx=v_ini*0.01/sqrt(1+(dydx(1,i))^2);
    X(1,i+1)=X(1,i)+dx;

    track(:,i)=[X(1,i);Y(1,i)];
end

direction=dydx;
end
% for i=2:N
%     disp(:,i)=[X(1,i)-X(1,i-1);Y(1,i)-Y(1,i-1)];
%     speed(i)=norm(disp(:,i))/0.01;                    
% end