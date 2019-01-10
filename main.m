close all
dt=0.01; tf=25; t=0:dt:tf; tfbar=30; t2bar=0:dt:tfbar; L=size(t2bar,2);
N=size(t,2); T=0.5; delta_T=0.002*T;
alpha=30;
Nbar=N;
tbar=t;
r_g=zeros(3,2,N);
for counter=1:3

    if counter==1
        V_x_ini=10;
    elseif counter==2
        V_x_ini=15; 
    else
        V_x_ini=19;
    end

udot=zeros(2,N);

% Car parameters 
lr=1.738; lf=1.105; m=2050; Iz=3344; Calpha=[-57500;-92500];      
accln_sat=50;
angle_sat=pi/2;

% Z is xdot,ydot,phidot, X, Y, phi
Z=zeros(6,N);
Z0=[25;0;0;0;0;1.13];
Z(:,1)=[V_x_ini;0;0;0;0;0];
%Z(:,1)=Z0;

%define reference r(1) is x_pos, r(2) is y_pos;

% V=zeros(2,L);M=N-1;
% 
% 
% for i=1:M/5
%     V(:,i)=[V_x_ini; 0];
% end
% 
% for i=1:M/10
%     V(:,i+M/5)=[sqrt(V_x_ini^2-(1.6*i/(M/10))^2);1.6*i/(M/10)];
% end
% 
% for i=1:M/10
%     V(:,i+3*M/10)=[sqrt(V_x_ini^2-(1.6-1.6*i/(M/10))^2);1.6-1.6*i/(M/10)];
% end
% 
% for i=1:3*M/25
%     V(:,i+2*M/5)=[sqrt(V_x_ini^2-(-2*i/(3*M/25))^2);-2*i/(3*M/25)];
% end    
% 
% for i=1:3*M/25
%     V(:,i+13*M/25)=[sqrt(V_x_ini^2-(-2+2*i/(3*M/25))^2);-2+2*i/(3*M/25)];
% end 
% 
% for i=16*M/25:M+1
%      V(:,i)=[V_x_ini;0];
% end 
% 
% for i=M+2:L
%     V(:,i)=[V_x_ini;0];
% end

r0=[0;0];
r=zeros(2,N);
%rextend=zeros(2,L);
r(:,1)=r0;

% for i=2:N
%     r(:,i)=r(:,i-1)+V(:,i)*dt;
% end    

[rextend, dydx]=trajectory(V_x_ini,tfbar);
%Define future reference, T time ahead
r=rextend(:,1:N);
for i=1:N
   ref(:,i)=rextend(:,i+T/dt); 
end    

% u is delta_f (front wheel angle wrt body) and a (acceleration)
u=zeros(2,N);
u(:,1)=[0;0];

% Track  vehicle motion
pos_tracker=zeros(2,N);

% derivative of z at all time instances
zdot=zeros(6,N);

g_u=zeros(2,N);
g_u_prime=zeros(2,2,N);
inv_matrix=zeros(2,2,N);


for i=2:N 
    if mod(i,1000)==0
        i
    end
    
    zdot(:,i)=dzdt(Z(:,i-1),u(:,i-1),lr,lf,m,Iz,Calpha);
  
%     for k=1:6  
%         if abs(zdot(k,i))>1
%             zdot(k,i)=zdot(k,i)/abs(zdot(k,i)); 
%         end
%     end    
     
    Z(:,i)=Z(:,i-1)+zdot(:,i)*dt;  
    pos_tracker(:,i)=[Z(4,i);Z(5,i)];
    
    %Future prediction
    [g_u(:,i),g_u_prime(:,:,i)]=g_rt(Z(:,i-1),u(:,i-1), T, delta_T, lr,lf,Calpha,m,Iz);
    
    if i==2
       g_u(:,1)=g_u(:,2); 
    end    
    
    inv_matrix(:,:,i)=inv(g_u_prime(:,:,i));
    
    u(:,i)=u(:,i-1)+alpha*(inv_matrix(:,:,i))*(ref(:,i)-g_u(:,i))*dt; 
    r_g(counter,:,i)=ref(:,i)-g_u(:,i);
    
    if i>N-T/(dt)
        %delta=u(:,N-T/(dt))/(T/(dt));
        %u(:,i)=u(:,i-1)-delta;
        u(:,i)=u(:,i-1);
    end    
    
    if abs(u(1,i))>angle_sat
        u(1,i)=u(1,i)*angle_sat/abs(u(1,i));
    end    
    
    if abs(u(2,i))>accln_sat
        u(2,i)=u(2,i)*accln_sat/abs(u(2,i));
    end
    
    udot(:,i-1)=(u(:,i)-u(:,i-1))/dt;
    
end    

lat_accl(counter,:)=u(2,:);
jerk(counter,:)=udot(2,:);


path_heading=mod(atan2(dydx(1,1:Nbar),1),2*pi);
actual_heading=mod(Z(6,1:Nbar),2*pi);

error_heading=mod(actual_heading-path_heading,2*pi);
for i=1:size(error_heading,2)
        if error_heading(i)>6
            error_heading(i)=error_heading(i)-2*pi;
        end 
end

head_error(counter,:)=180/pi*error_heading;


%[close_point,lateral_error_norm(counter,:),arc]=distance2curve(r',pos_tracker');
lateral_error_norm(counter,:)=abs(r(2,:)-pos_tracker(2,:));

error=pos_tracker-r;   
norm_error(counter,:)=vecnorm(error(:,1:Nbar));

position(counter,:,:)=pos_tracker(:,1:Nbar);
end

%% Plot of accln
figure (1)
plot(tbar,lat_accl(1,:),'LineWidth',1.5)
hold on
plot(tbar,lat_accl(2,:),'LineWidth',1.5)
plot(tbar,lat_accl(3,:),'LineWidth',1.5)

x1=xlabel('Time$~[s]$');
 y1=ylabel('$a_l~[m/s^2]$');
 set(x1,'Interpreter','latex')
 set(y1,'Interpreter','latex')
 
  leg1=legend('$10~m/s$','$15~m/s$','$19~m/s$');
 set(leg1,'Interpreter','latex')
 
 set(gcf, 'color', 'none');
set(gca, 'color', 'none');
hold off

pbaspect([2.5 1 1])
fig.PaperUnits = 'inches';
print('accln','-dsvg','-r0')
%% Plot of jerk
% figure (2)
% plot(tbar,jerk(1,:),'LineWidth',1.5)
% hold on
% plot(tbar,jerk(2,:),'LineWidth',1.5)
% plot(tbar,jerk(3,:),'LineWidth',1.5)
% %title('Longitudinal jerk vs time');
% 
% x1=xlabel('$time~[s]$');
%  y1=ylabel('$Jerk~[m/s^3]$');
%   set(x1,'Interpreter','latex')
%  set(y1,'Interpreter','latex')
%  
%  leg1=legend('$10~m/s$','$15~m/s$','$19~m/s$');
%  set(leg1,'Interpreter','latex')
%  set(gcf, 'color', 'none');
% set(gca, 'color', 'none');
% hold off
% 
% pbaspect([2.5 1 1])
% fig.PaperUnits = 'inches';
% print('jerk','-dsvg','-r0')
%% PLot of tracking error norm
% figure (3)
% 
% 
% plot (tbar,norm_error(1,:),'LineWidth',1.5);
% hold on
% plot (tbar,norm_error(2,:),'LineWidth',1.5);
% plot (tbar,norm_error(3,:),'LineWidth',1.5);
% %title('tracking error norm vs time');
% x1=xlabel('$time~[s]$');
%  y1=ylabel('Total $Error~[m]$');
%   set(x1,'Interpreter','latex')
%  set(y1,'Interpreter','latex')
%  
%  leg1=legend('$10~m/s$','$15~m/s$','$19~m/s$');
%  set(leg1,'Interpreter','latex')
%  
% set(gcf, 'color', 'none');
% set(gca, 'color', 'none');
% hold off
% 
% pbaspect([2.5 1 1])
% fig.PaperUnits = 'inches';
% print('total_error','-dsvg','-r0')
%% plot of heading error 
figure (4)

plot(tbar,head_error(1,:),'LineWidth',1.5)
hold on
plot(tbar,head_error(2,:),'LineWidth',1.5)
plot(tbar,head_error(3,:),'LineWidth',1.5)
%title('Heading error vs time')
x1=xlabel('Time$~[s]$');
y1=ylabel('$\Delta\psi~[^\circ]$');
 set(x1,'Interpreter','latex')
 set(y1,'Interpreter','latex')

leg1=legend('$10~m/s$','$15~m/s$','$19~m/s$');
 set(leg1,'Interpreter','latex')
     set(gcf, 'color', 'none');
set(gca, 'color', 'none');
hold off

pbaspect([2.5 1 1])
fig.PaperUnits = 'inches';
print('head_error','-dsvg','-r0')
%% PLot of lateral error
figure (5)

plot(tbar,lateral_error_norm(1,1:Nbar),'LineWidth',1.5)
hold on
plot(tbar,lateral_error_norm(2,1:Nbar),'LineWidth',1.5)
plot(tbar,lateral_error_norm(3,1:Nbar),'LineWidth',1.5)
%title ('Normal Error')
 x1=xlabel('Time$~[s]$');
 y1=ylabel('Lateral Error$~[m]$');
  set(x1,'Interpreter','latex')
 set(y1,'Interpreter','latex')
 
 leg1=legend('$10~m/s$','$15~m/s$','$19~m/s$');
 set(leg1,'Interpreter','latex')
     
set(gcf, 'color', 'none');
set(gca, 'color', 'none');
hold off

pbaspect([2.5 1 1])
fig.PaperUnits = 'inches';
print('lat_error','-dsvg','-r0')
 %% plot of path
 figure (6)

 plot (r(1,1:Nbar),r(2,1:Nbar),'LineWidth',1.5);
 hold on
 pos_15_x(1,:)=position(1,1,1:Nbar); pos_15_y(1,:)=position(1,2,1:Nbar);
 pos_25_x(1,:)=position(2,1,1:Nbar); pos_25_y(1,:)=position(2,2,1:Nbar);
 pos_35_x(1,:)=position(3,1,1:Nbar); pos_35_y(1,:)=position(3,2,1:Nbar);
 
 plot(pos_15_x,pos_15_y,'LineWidth',1.5);
 plot(pos_25_x,pos_25_y,'LineWidth',1.5);
 plot(pos_35_x,pos_35_y,'LineWidth',1.5);

x1=xlabel('$z_1~[m]$');
 y1=ylabel('$z_2~[m]$');
  set(x1,'Interpreter','latex')
 set(y1,'Interpreter','latex')
 
leg1=legend('Reference','$10~m/s$','$15~m/s$','$19~m/s$');
 set(leg1,'Interpreter','latex')
 
 set(gcf, 'color', 'none');
set(gca, 'color', 'none');
hold off

pbaspect([2.5 1 1])
fig.PaperUnits = 'inches';
print('path','-dsvg','-r0')
%% plot of r-g(u)
figure (8)
for i=1:N
    e1(1,i)=norm(r_g(1,:,i));
    e2(1,i)=norm(r_g(2,:,i));
    e3(1,i)=norm(r_g(3,:,i));
end

plot (t,e1,'LineWidth',1.5);
hold on
plot (t,e2,'LineWidth',1.5);
plot (t,e3,'LineWidth',1.5);

x1=xlabel('Time$~[s]$');
 y1=ylabel('Control Error$~[m]$');
  set(x1,'Interpreter','latex')
 set(y1,'Interpreter','latex')
 
leg1=legend('$10~m/s$','$15~m/s$','$19~m/s$');
 set(leg1,'Interpreter','latex')
 
 set(gcf, 'color', 'none');
set(gca, 'color', 'none');
hold off

pbaspect([2.5 1 1])
fig.PaperUnits = 'inches';
print('Control_error','-dsvg','-r0')

% Plot of accln    
% figure (1)
% plot(t,u(2,:))
% 
% title('longitudinal acceleration vs time');
% 
% legend('actual accleration');
% 
% 
% 
% % Plot of jerk
% figure (2)
% plot(t,udot(2,:))
% title('Longitudinal jerk vs time');
% 
% 
% % PLot of tracking error norm
% figure (3)
% error=pos_tracker-r;   
% 
% norm_error=vecnorm(error);
% 
% plot (t,norm_error);
% title('tracking error norm vs time');
% 
% 
% figure (4)
% 
% path_heading=mod(atan2(V(2,1:N),V(1,1:N)),2*pi);
% actual_heading=mod(Z(6,1:N),2*pi);
% 
% error_heading=mod(actual_heading-path_heading,2*pi);
% for i=1:size(error_heading,2)
%         if error_heading(i)>6
%             error_heading(i)=error_heading(i)-2*pi;
%         end 
% end
% 
% plot(t,error_heading)
% title('Heading error vs time');
% 
% 
% 
% figure (5)
% 
% [close_point,lateral_error_norm,arc]=distance2curve(r',pos_tracker');
% 
% plot(t,lateral_error_norm)
% title ('Normal Error')
%  
%     
% figure (6)
% plot(pos_tracker(1,:),pos_tracker(2,:), '--r');
% hold on
% plot (r(1,:),r(2,:),'b');
% title ('Tracking path')
% % xlim([-20 100]) 
% % ylim([-20 100])
% % 
% % s = plot(pos_tracker(1,1),pos_tracker(2,1),'o','MarkerFaceColor','black');      % bot 1
% % q = plot(r(1,1),r(2,1),'o','MarkerFaceColor','green');                          % ref
% % %r = 
% % 
% % for k = 2:N
% %     s.XData = pos_tracker(1,k);
% %     s.YData = pos_tracker(2,k);
% %       
% %     q.XData = r(1,k);
% %     q.YData = r(2,k);
% %     
% %     drawnow
% % end
%     
% 
% figure (7)
% plot(t,r(2,:))
% hold on
% plot(t,Z(5,:))
% legend('reference','actual');
% title ('Y position')



