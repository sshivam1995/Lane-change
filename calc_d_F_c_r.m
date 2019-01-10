function d_alpha_r=calc_d_F_c_r(z,u,C_alpha, lf,lr)

    d_alpha_r=zeros(8,1);
    
    d_alpha_r(1,1)= -1/(1+((z(2)-lr*z(3))/z(1))^2)*(-(z(2)-lr*z(3))/z(1)^2);
    d_alpha_r(2,1)= -1/(1+((z(2)-lr*z(3))/z(1))^2)*(1/z(1));
    d_alpha_r(3,1)= -1/(1+((z(2)-lr*z(3))/z(1))^2)*(-lr/z(1));
    d_alpha_r(4,1)= 0;
    d_alpha_r(5,1)= 0;
    d_alpha_r(6,1)= 0;
    d_alpha_r(7,1)= 0;
    d_alpha_r(8,1)= 0;
    
    d_alpha_r=d_alpha_r*(-C_alpha(2));
end