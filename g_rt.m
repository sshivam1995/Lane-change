function  [g_u,g_u_prime]= g_rt(z0,u,T,delta_t,lr,lf,C_alpha,m,Iz)

C=[0 0 0 1 0 0; 0 0 0 0 1 0];

M=fix(T/delta_t);

z=zeros(6,M);
z(:,1)=z0;

dzdu=zeros(6,2,M);

for i=2:M
    zdot(:,i)=dzdt(z(:,i-1),u,lr,lf,m,Iz,C_alpha);
    z(:,i)=z(:,i-1)+zdot(:,i)*delta_t;
    
    dfdz=calcdfdz(z(:,i),u,lr,lf,C_alpha,m,Iz);
    dfdu=calcdfdu(z(:,i),u,lr,lf,C_alpha,m,Iz);
    
    dzdu(:,:,i)=dzdu(:,:,i-1)+(dfdz*dzdu(:,:,i-1)+dfdu)*delta_t;     
end

g_u=C*z(:,M);
g_u_prime=C*dzdu(:,:,M);

end




