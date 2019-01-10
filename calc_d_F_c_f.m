function d_alpha_f=calc_d_F_c_f(z,u,C_alpha, lf, lr)
    d_alpha_f=zeros(8,1);
    
    d_alpha_f(1,1)= -1/(1+((z(2)+lf*z(3))/z(1))^2)*(-(z(2)+lf*z(3))/z(1)^2);
    d_alpha_f(2,1)= -1/(1+((z(2)+lf*z(3))/z(1))^2)*(1/z(1));
    d_alpha_f(3,1)= -1/(1+((z(2)+lf*z(3))/z(1))^2)*(lf/z(1));
    d_alpha_f(4,1)= 0;
    d_alpha_f(5,1)= 0;
    d_alpha_f(6,1)= 0;
    d_alpha_f(7,1)= 1;
    d_alpha_f(8,1)= 0;
    
    d_alpha_f=d_alpha_f*(-C_alpha(1));
end