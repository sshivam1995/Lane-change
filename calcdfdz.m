function dfdz=calcdfdz(z,u,lr,lf,C_alpha,m,Iz)

    d_F_c_f=calc_d_F_c_f(z,u,C_alpha,lf,lr);
    d_F_c_r=calc_d_F_c_r(z,u,C_alpha,lf,lr);
    
    dfdz=zeros(6,6);
    
    dfdz(1,1)= 0; 
    dfdz(1,2)= z(3); 
    dfdz(1,3)= z(2); 
    dfdz(1,4)= 0;
    dfdz(1,5)= 0;
    dfdz(1,6)= 0;
    
    dfdz(2,1)= -z(3)+(2/m)*(d_F_c_f(1)*cos(u(1))+d_F_c_r(1));
    dfdz(2,2)= (2/m)*(d_F_c_f(2)*cos(u(1))+d_F_c_r(2));
    dfdz(2,3)= -z(1)+(2/m)*(d_F_c_f(3)*cos(u(1))+d_F_c_r(3));
    dfdz(2,4)= (2/m)*(d_F_c_f(4)*cos(u(1))+d_F_c_r(4));
    dfdz(2,5)= (2/m)*(d_F_c_f(5)*cos(u(1))+d_F_c_r(5));
    dfdz(2,6)= (2/m)*(d_F_c_f(6)*cos(u(1))+d_F_c_r(6));
    
    dfdz(3,1)= (2/Iz)*(lf*d_F_c_f(1)-lr*d_F_c_r(1));
    dfdz(3,2)= (2/Iz)*(lf*d_F_c_f(2)-lr*d_F_c_r(2));
    dfdz(3,3)= (2/Iz)*(lf*d_F_c_f(3)-lr*d_F_c_r(3));
    dfdz(3,4)= (2/Iz)*(lf*d_F_c_f(4)-lr*d_F_c_r(4));
    dfdz(3,5)= (2/Iz)*(lf*d_F_c_f(5)-lr*d_F_c_r(5));
    dfdz(3,6)= (2/Iz)*(lf*d_F_c_f(6)-lr*d_F_c_r(6));
    
    dfdz(4,1)= cos(z(6));
    dfdz(4,2)= -sin(z(6));
    dfdz(4,3)= 0;
    dfdz(4,4)= 0;
    dfdz(4,5)= 0;
    dfdz(4,6)= -z(1)*sin(z(6))-z(2)*cos(z(6));
    
    dfdz(5,1)= sin(z(6));
    dfdz(5,2)= cos(z(6));
    dfdz(5,3)= 0;
    dfdz(5,4)= 0;
    dfdz(5,5)= 0;
    dfdz(5,6)= z(1)*cos(z(6))-z(2)*sin(z(6));
    
    dfdz(6,1)= 0;
    dfdz(6,2)= 0;
    dfdz(6,3)= 1;
    dfdz(6,4)= 0;
    dfdz(6,5)= 0;
    dfdz(6,6)= 0;
    
end    