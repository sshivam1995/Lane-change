function dfdu=calcdfdu(z,u,lr,lf,C_alpha,m,Iz)

    alpha_f = u(1)-atan((z(2)+lf*z(3))/z(1));        % delta-theta_v_f
    alpha_r = -atan((z(2)-lr*z(3))/z(1));           % -theta_v_r
    
    F_c_f=-C_alpha(1)*alpha_f;
    F_c_r=-C_alpha(2)*alpha_r;

    d_F_c_f=calc_d_F_c_f(z,u,C_alpha,lf,lr);
    d_F_c_r=calc_d_F_c_r(z,u,C_alpha,lf,lr);
    
    dfdu=zeros(6,2);
    
    dfdu(1,1)= 0; 
    dfdu(1,2)= 1; 
    
    dfdu(2,1)= (2/m)*(-F_c_f*sin(u(1))+d_F_c_f(7)*cos(u(1))+d_F_c_r(7)); 
    dfdu(2,2)= (2/m)*(d_F_c_f(8)*cos(u(1))+d_F_c_r(8));
    
    dfdu(3,1)= (2/Iz)*(lf*d_F_c_f(7)-lr*d_F_c_r(7));
    dfdu(3,2)= (2/Iz)*(lf*d_F_c_f(8)-lr*d_F_c_r(8));
    
    dfdu(4,1)= 0;
    dfdu(4,2)= 0;
    
    dfdu(5,1)= 0;
    dfdu(5,2)= 0;
    
    dfdu(6,1)= 0;
    dfdu(6,2)= 0;
    
end