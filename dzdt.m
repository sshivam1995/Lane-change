function zdot=dzdt(z,u,lr,lf,m,Iz,C_alpha)
    
    alpha_f = u(1)-atan((z(2)+lf*z(3))/z(1));        % delta-theta_v_f
    alpha_r = -atan((z(2)-lr*z(3))/z(1));           % -theta_v_r
    
    F_c_f=-C_alpha(1)*alpha_f;
    F_c_r=-C_alpha(2)*alpha_r;
    
    %beta=atan(lr/(lf+lr)*tan(u(1)));
    
    zdot(1,1)= z(3)*z(2)+u(2);
    zdot(2,1)= -z(3)*z(1)+(2/m)*(F_c_f*cos(u(1))+F_c_r);
    zdot(3,1)= (2/Iz)*(lf*F_c_f*cos(z(6))-lr*F_c_r);
    zdot(4,1)= z(1)*cos(z(6))-z(2)*sin(z(6));
    zdot(5,1)= z(1)*sin(z(6))+z(2)*cos(z(6));
    zdot(6,1)= z(3);
    
end