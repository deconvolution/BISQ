function [ux,uz,p,uxt,uzt,pt,recx,recz,recp]=poroelastic_solver_PLM(dims,C,Ks,rhos,Kf,rhof,eta,phi,k1,k3,R1,R3,rhoa,ome,lp,d0,vb)
%% empty matrices
ux=zeros(dims.nz,dims.nx,length(ome));
uz=zeros(dims.nz,dims.nx,length(ome));
p=zeros(dims.nz,dims.nx,length(ome));
s=zeros(dims.nz*dims.nx*3,1);
%% empty AA
AA=sparse(dims.nz*dims.nx*3,dims.nz*dims.nx*3);
%%
tic;
for l=2:length(ome)
    %% Frequency dependent parameters
    alp1=1-(C(1,1)+C(1,2)+C(1,3))./3./Ks;
    alp3=1-(2*C(1,3)+C(3,3))/3./Ks;
    
    % problem
    bet=phi./Kf+(1-phi)./Ks-(2.*(C(1,1)+C(1,2)+2*C(1,3))+C(3,3))/9./Ks.^2;
    
    K1ome=1i.*phi./ome(l)./rhof.*((rhoa./rhof+phi)./phi+1i.*eta.*phi./k1./rhof./ome(l)).^(-1);
    K3ome=1i.*phi./ome(l)./rhof.*((rhoa./rhof+phi)./phi+1i.*eta.*phi./k3./rhof./ome(l)).^(-1);
    
    gam1=sqrt(rhof.*ome(l).^2./(phi./bet).*((rhoa./rhof+phi)./phi+1i.*eta.*phi./k1./rhof./ome(l)));
    gam3=sqrt(rhof.*ome(l).^2./(phi./bet).*((rhoa./rhof+phi)./phi+1i.*eta.*phi./k3./rhof./ome(l)));
    gam1(isnan(gam1))=0;
    gam3(isnan(gam3))=0;
    
    s1=1-2*besselj(1,gam1.*R1)./gam1./R1./besselj(0,gam1.*R1);
    s3=1-2*besselj(1,gam3.*R3)./gam3./R3./besselj(0,gam3.*R3);
    s1(isnan(s1))=0;
    s3(isnan(s3))=0;
    
    rho=rhos.*(1-phi)+phi.*rhof;
    
    the1=-K1ome./1i./ome(l);
    the3=-K3ome./1i./ome(l);
    
    barrho1=rho+ome(l).^2.*rhof.^2.*the1;
    barrho3=rho+ome(l).^2.*rhof.^2.*the3;
    
    baralp1=alp1+ome(l).^2.*rhof.*the1;
    baralp3=alp3+ome(l).^2.*rhof.*the3;
    %% assign d
    dx=sparse(dims.nz,dims.nx);
    dz=sparse(dims.nz,dims.nx);
    % top
    for i=2:lp+1
        for j=lp+2:dims.nx-lp-1
            dx(i,j)=0;
            dz(i,j)=d0*((lp+2-i)/lp)^2;
        end
    end
    % bottom
    for i=dims.nz-lp:dims.nz-1
        for j=lp+2:dims.nx-lp-1
            dx(i,j)=0;
            dz(i,j)=d0*((i+1-(dims.nz-lp))/lp)^2;
        end
    end
    % left
    for i=lp+2:dims.nz-lp-1
        for j=2:lp+1
            dz(i,j)=0;
            dx(i,j)=d0*((lp+2-j)/lp)^2;
        end
    end
    % right
    for i=lp+1:dims.nz-lp-1
        for j=dims.nx-lp:dims.nx-1
            dz(i,j)=0;
            dx(i,j)=d0*((j+1-(dims.nx-lp))/lp)^2;
        end
    end
    % upper left
    for i=2:lp+1
        for j=2:lp+1
            dz(i,j)=d0*((lp+2-i)/lp)^2;
            dx(i,j)=d0*((lp+2-j)/lp)^2;
        end
    end
    % upper right
    for i=2:lp+1
        for j=dims.nx-lp:dims.nx-1
            dz(i,j)=d0*((lp+2-i)/lp)^2;
            dx(i,j)=d0*((j+1-(dims.nx-lp))/lp)^2;
        end
    end
    % lower left
    for i=dims.nz-lp:dims.nz-1
        for j=2:lp+1
            dz(i,j)=d0*((i+1-(dims.nz-lp))/lp)^2;
            dx(i,j)=d0*((lp+2-j)/lp)^2;
        end
    end
    % lower right
    for i=dims.nz-lp:dims.nz-1
        for j=dims.nx-lp:dims.nx-1
            dz(i,j)=d0*((i+1-(dims.nz-lp))/lp)^2;
            dx(i,j)=d0*((j+1-(dims.nx-lp))/lp)^2;
        end
    end
    gx=1+dx/1i./ome(l);
    gz=1+dz/1i./ome(l);
    gx_x=diff(gx,1,2)/dims.dh;
    gz_z=diff(gz,1,1)/dims.dh;
    gx_x=[-gx_x(:,1:lp+1),zeros([size(gx_x,1),1]),gx_x(:,lp+2:end)];
    gx_x(:,1)=0;
    gx_x(:,end)=0;
    gz_z=[-gz_z(1:lp+1,:);zeros([1,size(gz_z,2)]);gz_z(lp+2:end,:)];
    gz_z(1,:)=0;
    gz_z(end,:)=0;
    %% domain interior
    for i=lp+2:dims.nz-1
        for j=lp+2:dims.nx-1
            %% Eq. 1
            % ux
            % i,j
            AA(((j-1)*dims.nz+i-1)*3+1,((j-1)*dims.nz+i-1)*3+1)=ome(l)^2*barrho1(i,j)-2*C(1,1,i,j)/gx(i,j)^2/dims.dh^2-2*C(4,4,i,j)/gz(i,j)^2/dims.dh^2;
            % i,j+1
            AA(((j-1)*dims.nz+i-1)*3+1,((j)*dims.nz+i-1)*3+1)=-C(1,1,i,j)/2/gx(i,j)^3/dims.dh*gx_x(i,j)+C(1,1,i,j)/gx(i,j)^2/dims.dh^2;
            % i,j-1
            AA(((j-1)*dims.nz+i-1)*3+1,((j-2)*dims.nz+i-1)*3+1)=C(1,1,i,j)/2/gx(i,j)^3/dims.dh*gx_x(i,j)+C(1,1,i,j)/gx(i,j)^2/dims.dh^2;
            % i+1,j
            AA(((j-1)*dims.nz+i-1)*3+1,((j-1)*dims.nz+i)*3+1)=-C(4,4,i,j)/2/gz(i,j)^3/dims.dh*gz_z(i,j)+C(4,4,i,j)/gz(i,j)^2/dims.dh^2;
            % i-1,j
            AA(((j-1)*dims.nz+i-1)*3+1,((j-1)*dims.nz+i-2)*3+1)=C(4,4,i,j)/2/gz(i,j)^3/dims.dh*gz_z(i,j)+C(4,4,i,j)/gz(i,j)^2/dims.dh^2;
            % uz
            % i+1,j+1
            AA(((j-1)*dims.nz+i-1)*3+1,((j)*dims.nz+i)*3+1+1)=(C(1,3,i,j)+C(4,4,i,j))/4/gx(i,j)/gz(i,j)/dims.dh^2;
            % i-1,j+1
            AA(((j-1)*dims.nz+i-1)*3+1,((j)*dims.nz+i-2)*3+1+1)=-(C(1,3,i,j)+C(4,4,i,j))/4/gx(i,j)/gz(i,j)/dims.dh^2;
            % i+1,j-1
            AA(((j-1)*dims.nz+i-1)*3+1,((j-2)*dims.nz+i)*3+1+1)=-(C(1,3,i,j)+C(4,4,i,j))/4/gx(i,j)/gz(i,j)/dims.dh^2;
            % i-1,j-1
            AA(((j-1)*dims.nz+i-1)*3+1,((j-2)*dims.nz+i-2)*3+1+1)=(C(1,3,i,j)+C(4,4,i,j))/4/gx(i,j)/gz(i,j)/dims.dh^2;
            % p
            % i,j+1
            AA(((j-1)*dims.nz+i-1)*3+1,((j)*dims.nz+i-1)*3+1+2)=-baralp1(i,j)/2/gx(i,j)/dims.dh;
            % i,j-1
            AA(((j-1)*dims.nz+i-1)*3+1,((j-2)*dims.nz+i-1)*3+1+2)=baralp1(i,j)/2/gx(i,j)/dims.dh;
            
            %% Eq. 2
            % ux
            % i+1,j+1
            AA(((j-1)*dims.nz+i-1)*3+1+1,((j)*dims.nz+i)*3+1)=(C(1,3,i,j)+C(4,4,i,j))/4/gx(i,j)/gz(i,j)/dims.dh^2;
            % i-1,j+1
            AA(((j-1)*dims.nz+i-1)*3+1+1,((j)*dims.nz+i-2)*3+1)=-(C(1,3,i,j)+C(4,4,i,j))/4/gx(i,j)/gz(i,j)/dims.dh^2;
            % i+1,j-1
            AA(((j-1)*dims.nz+i-1)*3+1+1,((j-2)*dims.nz+i)*3+1)=-(C(1,3,i,j)+C(4,4,i,j))/4/gx(i,j)/gz(i,j)/dims.dh^2;
            % i-1,j-1
            AA(((j-1)*dims.nz+i-1)*3+1+1,((j-2)*dims.nz+i-2)*3+1)=(C(1,3,i,j)+C(4,4,i,j))/4/gx(i,j)/gz(i,j)/dims.dh^2;
            % uz
            % i,j
            AA(((j-1)*dims.nz+i-1)*3+1+1,((j-1)*dims.nz+i-1)*3+1+1)=ome(l)^2*barrho3(i,j)-2*C(4,4,i,j)/gx(i,j)^2/dims.dh^2-2*C(3,3,i,j)/gz(i,j)^2/dims.dh^2;
            % i,j+1
            AA(((j-1)*dims.nz+i-1)*3+1+1,((j)*dims.nz+i-1)*3+1+1)=-C(4,4,i,j)/2/gx(i,j)^3/dims.dh*gx_x(i,j)+C(4,4,i,j)/gx(i,j)^2/dims.dh^2;
            % i,j-1
            AA(((j-1)*dims.nz+i-1)*3+1+1,((j-2)*dims.nz+i-1)*3+1+1)=C(4,4,i,j)/2/gx(i,j)^3/dims.dh*gx_x(i,j)+C(4,4,i,j)/gx(i,j)^2/dims.dh^2;
            % i+1,j
            AA(((j-1)*dims.nz+i-1)*3+1+1,((j-1)*dims.nz+i)*3+1+1)=-C(3,3,i,j)/2/gz(i,j)^3/dims.dh*gz_z(i,j)+C(3,3,i,j)/gz(i,j)^2/dims.dh^2;
            % i-1,j
            AA(((j-1)*dims.nz+i-1)*3+1+1,((j-1)*dims.nz+i-2)*3+1+1)=C(3,3,i,j)/2/gz(i,j)^3/dims.dh*gz_z(i,j)+C(3,3,i,j)/gz(i,j)^2/dims.dh^2;
            % p
            % i+1,j
            AA(((j-1)*dims.nz+i-1)*3+1+1,((j-1)*dims.nz+i)*3+1+2)=-baralp3(i,j)/2/gz(i,j)/dims.dh;
            % i-1,j
            AA(((j-1)*dims.nz+i-1)*3+1+1,((j-1)*dims.nz+i-2)*3+1+2)=baralp3(i,j)/2/gz(i,j)/dims.dh;
            %% Eq. 3
            % ux
            % i,j+1
            AA(((j-1)*dims.nz+i-1)*3+1+2,((j)*dims.nz+i-1)*3+1)=-s3(i,j)*baralp1(i,j)/2/gx(i,j)/dims.dh;
            % i,j-1
            AA(((j-1)*dims.nz+i-1)*3+1+2,((j-2)*dims.nz+i-1)*3+1)=s3(i,j)*baralp1(i,j)/2/gx(i,j)/dims.dh;
            % uz
            % i+1,j
            AA(((j-1)*dims.nz+i-1)*3+1+2,((j-1)*dims.nz+i)*3+1+1)=-s1(i,j)*baralp3(i,j)/2/gz(i,j)/dims.dh;
            % i-1,j
            AA(((j-1)*dims.nz+i-1)*3+1+2,((j-1)*dims.nz+i-2)*3+1+1)=s1(i,j)*baralp3(i,j)/2/gz(i,j)/dims.dh;
            % p
            % i,j
            AA(((j-1)*dims.nz+i-1)*3+1+2,((j-1)*dims.nz+i-1)*3+1+2)=-bet(i,j)-2*s1(i,j)*the3(i,j)/gz(i,j)^2/dims.dh^2-2*s3(i,j)*the1(i,j)/gx(i,j)^2/dims.dh^2;
            % i,j+1
            AA(((j-1)*dims.nz+i-1)*3+1+2,((j)*dims.nz+i-1)*3+1+2)=-s3(i,j)*the1(i,j)/2/gx(i,j)^3/dims.dh*gx_x(i,j)+s3(i,j)*the1(i,j)/gx(i,j)^2/dims.dh^2;
            % i,j-1
            AA(((j-1)*dims.nz+i-1)*3+1+2,((j-2)*dims.nz+i-1)*3+1+2)=s3(i,j)*the1(i,j)/2/gx(i,j)^3/dims.dh*gx_x(i,j)+s3(i,j)*the1(i,j)/gx(i,j)^2/dims.dh^2;
            % i+1,j
            AA(((j-1)*dims.nz+i-1)*3+1+2,((j-1)*dims.nz+i)*3+1+2)=-s1(i,j)*the3(i,j)/2/gz(i,j)^3/dims.dh*gz_z(i,j)+s1(i,j)*the3(i,j)/gz(i,j)^2/dims.dh^2;
            % i-1,j
            AA(((j-1)*dims.nz+i-1)*3+1+2,((j-1)*dims.nz+i-2)*3+1+2)=s1(i,j)*the3(i,j)/2/gz(i,j)^3/dims.dh*gz_z(i,j)+s1(i,j)*the3(i,j)/gz(i,j)^2/dims.dh^2;
        end
    end
    %% top local
    for i=2:lp+1
        for j=lp+2:dims.nx-lp-1
            %% Eq. 1
            % ux
            % i,j
            AA(((j-1)*dims.nz+i-1)*3+1,((j-1)*dims.nz+i-1)*3+1)=ome(l)^2*barrho1(i,j)-2*C(1,1,i,j)/gx(i,j)^2/dims.dh^2-2*C(4,4,i,j)/gz(i,j)^2/dims.dh^2;
            % i,j+1
            AA(((j-1)*dims.nz+i-1)*3+1,((j)*dims.nz+i-1)*3+1)=-C(1,1,i,j)/2/gx(i,j)^3/dims.dh*gx_x(i,j)+C(1,1,i,j)/gx(i,j)^2/dims.dh^2;
            % i,j-1
            AA(((j-1)*dims.nz+i-1)*3+1,((j-2)*dims.nz+i-1)*3+1)=C(1,1,i,j)/2/gx(i,j)^3/dims.dh*gx_x(i,j)+C(1,1,i,j)/gx(i,j)^2/dims.dh^2;
            % i+1,j
            AA(((j-1)*dims.nz+i-1)*3+1,((j-1)*dims.nz+i)*3+1)=C(4,4,i,j)/2/gz(i,j)^3/dims.dh*gz_z(i,j)+C(4,4,i,j)/gz(i,j)^2/dims.dh^2;
            % i-1,j
            AA(((j-1)*dims.nz+i-1)*3+1,((j-1)*dims.nz+i-2)*3+1)=-C(4,4,i,j)/2/gz(i,j)^3/dims.dh*gz_z(i,j)+C(4,4,i,j)/gz(i,j)^2/dims.dh^2;
            % uz
            % i+1,j+1
            AA(((j-1)*dims.nz+i-1)*3+1,((j)*dims.nz+i)*3+1+1)=(C(1,3,i,j)+C(4,4,i,j))/4/gx(i,j)/gz(i,j)/dims.dh^2;
            % i-1,j+1
            AA(((j-1)*dims.nz+i-1)*3+1,((j)*dims.nz+i-2)*3+1+1)=-(C(1,3,i,j)+C(4,4,i,j))/4/gx(i,j)/gz(i,j)/dims.dh^2;
            % i+1,j-1
            AA(((j-1)*dims.nz+i-1)*3+1,((j-2)*dims.nz+i)*3+1+1)=-(C(1,3,i,j)+C(4,4,i,j))/4/gx(i,j)/gz(i,j)/dims.dh^2;
            % i-1,j-1
            AA(((j-1)*dims.nz+i-1)*3+1,((j-2)*dims.nz+i-2)*3+1+1)=(C(1,3,i,j)+C(4,4,i,j))/4/gx(i,j)/gz(i,j)/dims.dh^2;
            % p
            % i,j+1
            AA(((j-1)*dims.nz+i-1)*3+1,((j)*dims.nz+i-1)*3+1+2)=-baralp1(i,j)/2/gx(i,j)/dims.dh;
            % i,j-1
            AA(((j-1)*dims.nz+i-1)*3+1,((j-2)*dims.nz+i-1)*3+1+2)=baralp1(i,j)/2/gx(i,j)/dims.dh;
            %% Eq. 2
            % ux
            % i+1,j+1
            AA(((j-1)*dims.nz+i-1)*3+1+1,((j)*dims.nz+i)*3+1)=(C(1,3,i,j)+C(4,4,i,j))/4/gx(i,j)/gz(i,j)/dims.dh^2;
            % i-1,j+1
            AA(((j-1)*dims.nz+i-1)*3+1+1,((j)*dims.nz+i-2)*3+1)=-(C(1,3,i,j)+C(4,4,i,j))/4/gx(i,j)/gz(i,j)/dims.dh^2;
            % i+1,j-1
            AA(((j-1)*dims.nz+i-1)*3+1+1,((j-2)*dims.nz+i)*3+1)=-(C(1,3,i,j)+C(4,4,i,j))/4/gx(i,j)/gz(i,j)/dims.dh^2;
            % i-1,j-1
            AA(((j-1)*dims.nz+i-1)*3+1+1,((j-2)*dims.nz+i-2)*3+1)=(C(1,3,i,j)+C(4,4,i,j))/4/gx(i,j)/gz(i,j)/dims.dh^2;
            % uz
            % i,j
            AA(((j-1)*dims.nz+i-1)*3+1+1,((j-1)*dims.nz+i-1)*3+1+1)=ome(l)^2*barrho3(i,j)-2*C(4,4,i,j)/gx(i,j)^2/dims.dh^2-2*C(3,3,i,j)/gz(i,j)^2/dims.dh^2;
            % i,j+1
            AA(((j-1)*dims.nz+i-1)*3+1+1,((j)*dims.nz+i-1)*3+1+1)=-C(4,4,i,j)/2/gx(i,j)^3/dims.dh*gx_x(i,j)+C(4,4,i,j)/gx(i,j)^2/dims.dh^2;
            % i,j-1
            AA(((j-1)*dims.nz+i-1)*3+1+1,((j-2)*dims.nz+i-1)*3+1+1)=C(4,4,i,j)/2/gx(i,j)^3/dims.dh*gx_x(i,j)+C(4,4,i,j)/gx(i,j)^2/dims.dh^2;
            % i+1,j
            AA(((j-1)*dims.nz+i-1)*3+1+1,((j-1)*dims.nz+i)*3+1+1)=C(3,3,i,j)/2/gz(i,j)^3/dims.dh*gz_z(i,j)+C(3,3,i,j)/gz(i,j)^2/dims.dh^2;
            % i-1,j
            AA(((j-1)*dims.nz+i-1)*3+1+1,((j-1)*dims.nz+i-2)*3+1+1)=-C(3,3,i,j)/2/gz(i,j)^3/dims.dh*gz_z(i,j)+C(3,3,i,j)/gz(i,j)^2/dims.dh^2;
            % p
            % i+1,j
            AA(((j-1)*dims.nz+i-1)*3+1+1,((j-1)*dims.nz+i)*3+1+2)=-baralp3(i,j)/2/gz(i,j)/dims.dh;
            % i-1,j
            AA(((j-1)*dims.nz+i-1)*3+1+1,((j-1)*dims.nz+i-2)*3+1+2)=baralp3(i,j)/2/gz(i,j)/dims.dh;
            %% Eq. 3
            % ux
            % i,j+1
            AA(((j-1)*dims.nz+i-1)*3+1+2,((j)*dims.nz+i-1)*3+1)=-s3(i,j)*baralp1(i,j)/2/gx(i,j)/dims.dh;
            % i,j-1
            AA(((j-1)*dims.nz+i-1)*3+1+2,((j-2)*dims.nz+i-1)*3+1)=s3(i,j)*baralp1(i,j)/2/gx(i,j)/dims.dh;
            % uz
            % i+1,j
            AA(((j-1)*dims.nz+i-1)*3+1+2,((j-1)*dims.nz+i)*3+1+1)=-s1(i,j)*baralp3(i,j)/2/gz(i,j)/dims.dh;
            % i-1,j
            AA(((j-1)*dims.nz+i-1)*3+1+2,((j-1)*dims.nz+i-2)*3+1+1)=s1(i,j)*baralp3(i,j)/2/gz(i,j)/dims.dh;
            % p
            % i,j
            AA(((j-1)*dims.nz+i-1)*3+1+2,((j-1)*dims.nz+i-1)*3+1+2)=-bet(i,j)-2*s1(i,j)*the3(i,j)/gz(i,j)^2/dims.dh^2-2*s3(i,j)*the1(i,j)/gx(i,j)^2/dims.dh^2;
            % i,j+1
            AA(((j-1)*dims.nz+i-1)*3+1+2,((j)*dims.nz+i-1)*3+1+2)=-s3(i,j)*the1(i,j)/2/gx(i,j)^3/dims.dh*gx_x(i,j)+s3(i,j)*the1(i,j)/gx(i,j)^2/dims.dh^2;
            % i,j-1
            AA(((j-1)*dims.nz+i-1)*3+1+2,((j-2)*dims.nz+i-1)*3+1+2)=s3(i,j)*the1(i,j)/2/gx(i,j)^3/dims.dh*gx_x(i,j)+s3(i,j)*the1(i,j)/gx(i,j)^2/dims.dh^2;
            % i+1,j
            AA(((j-1)*dims.nz+i-1)*3+1+2,((j-1)*dims.nz+i)*3+1+2)=s1(i,j)*the3(i,j)/2/gz(i,j)^3/dims.dh*gz_z(i,j)+s1(i,j)*the3(i,j)/gz(i,j)^2/dims.dh^2;
            % i-1,j
            AA(((j-1)*dims.nz+i-1)*3+1+2,((j-1)*dims.nz+i-2)*3+1+2)=-s1(i,j)*the3(i,j)/2/gz(i,j)^3/dims.dh*gz_z(i,j)+s1(i,j)*the3(i,j)/gz(i,j)^2/dims.dh^2;
        end
    end
    %% left local
    for i=lp+2:dims.nz-lp-1
        for j=2:lp+1
            %% Eq. 1
            % ux
            % i,j
            AA(((j-1)*dims.nz+i-1)*3+1,((j-1)*dims.nz+i-1)*3+1)=ome(l)^2*barrho1(i,j)-2*C(1,1,i,j)/gx(i,j)^2/dims.dh^2-2*C(4,4,i,j)/gz(i,j)^2/dims.dh^2;
            % i,j+1
            AA(((j-1)*dims.nz+i-1)*3+1,((j)*dims.nz+i-1)*3+1)=C(1,1,i,j)/2/gx(i,j)^3/dims.dh*gx_x(i,j)+C(1,1,i,j)/gx(i,j)^2/dims.dh^2;
            % i,j-1
            AA(((j-1)*dims.nz+i-1)*3+1,((j-2)*dims.nz+i-1)*3+1)=-C(1,1,i,j)/2/gx(i,j)^3/dims.dh*gx_x(i,j)+C(1,1,i,j)/gx(i,j)^2/dims.dh^2;
            % i+1,j
            AA(((j-1)*dims.nz+i-1)*3+1,((j-1)*dims.nz+i)*3+1)=-C(4,4,i,j)/2/gz(i,j)^3/dims.dh*gz_z(i,j)+C(4,4,i,j)/gz(i,j)^2/dims.dh^2;
            % i-1,j
            AA(((j-1)*dims.nz+i-1)*3+1,((j-1)*dims.nz+i-2)*3+1)=C(4,4,i,j)/2/gz(i,j)^3/dims.dh*gz_z(i,j)+C(4,4,i,j)/gz(i,j)^2/dims.dh^2;
            % uz
            % i+1,j+1
            AA(((j-1)*dims.nz+i-1)*3+1,((j)*dims.nz+i)*3+1+1)=(C(1,3,i,j)+C(4,4,i,j))/4/gx(i,j)/gz(i,j)/dims.dh^2;
            % i-1,j+1
            AA(((j-1)*dims.nz+i-1)*3+1,((j)*dims.nz+i-2)*3+1+1)=-(C(1,3,i,j)+C(4,4,i,j))/4/gx(i,j)/gz(i,j)/dims.dh^2;
            % i+1,j-1
            AA(((j-1)*dims.nz+i-1)*3+1,((j-2)*dims.nz+i)*3+1+1)=-(C(1,3,i,j)+C(4,4,i,j))/4/gx(i,j)/gz(i,j)/dims.dh^2;
            % i-1,j-1
            AA(((j-1)*dims.nz+i-1)*3+1,((j-2)*dims.nz+i-2)*3+1+1)=(C(1,3,i,j)+C(4,4,i,j))/4/gx(i,j)/gz(i,j)/dims.dh^2;
            % p
            % i,j+1
            AA(((j-1)*dims.nz+i-1)*3+1,((j)*dims.nz+i-1)*3+1+2)=-baralp1(i,j)/2/gx(i,j)/dims.dh;
            % i,j-1
            AA(((j-1)*dims.nz+i-1)*3+1,((j-2)*dims.nz+i-1)*3+1+2)=baralp1(i,j)/2/gx(i,j)/dims.dh;
            
            %% Eq. 2
            % ux
            % i+1,j+1
            AA(((j-1)*dims.nz+i-1)*3+1+1,((j)*dims.nz+i)*3+1)=(C(1,3,i,j)+C(4,4,i,j))/4/gx(i,j)/gz(i,j)/dims.dh^2;
            % i-1,j+1
            AA(((j-1)*dims.nz+i-1)*3+1+1,((j)*dims.nz+i-2)*3+1)=-(C(1,3,i,j)+C(4,4,i,j))/4/gx(i,j)/gz(i,j)/dims.dh^2;
            % i+1,j-1
            AA(((j-1)*dims.nz+i-1)*3+1+1,((j-2)*dims.nz+i)*3+1)=-(C(1,3,i,j)+C(4,4,i,j))/4/gx(i,j)/gz(i,j)/dims.dh^2;
            % i-1,j-1
            AA(((j-1)*dims.nz+i-1)*3+1+1,((j-2)*dims.nz+i-2)*3+1)=(C(1,3,i,j)+C(4,4,i,j))/4/gx(i,j)/gz(i,j)/dims.dh^2;
            % uz
            % i,j
            AA(((j-1)*dims.nz+i-1)*3+1+1,((j-1)*dims.nz+i-1)*3+1+1)=ome(l)^2*barrho3(i,j)-2*C(4,4,i,j)/gx(i,j)^2/dims.dh^2-2*C(3,3,i,j)/gz(i,j)^2/dims.dh^2;
            % i,j+1
            AA(((j-1)*dims.nz+i-1)*3+1+1,((j)*dims.nz+i-1)*3+1+1)=C(4,4,i,j)/2/gx(i,j)^3/dims.dh*gx_x(i,j)+C(4,4,i,j)/gx(i,j)^2/dims.dh^2;
            % i,j-1
            AA(((j-1)*dims.nz+i-1)*3+1+1,((j-2)*dims.nz+i-1)*3+1+1)=-C(4,4,i,j)/2/gx(i,j)^3/dims.dh*gx_x(i,j)+C(4,4,i,j)/gx(i,j)^2/dims.dh^2;
            % i+1,j
            AA(((j-1)*dims.nz+i-1)*3+1+1,((j-1)*dims.nz+i)*3+1+1)=-C(3,3,i,j)/2/gz(i,j)^3/dims.dh*gz_z(i,j)+C(3,3,i,j)/gz(i,j)^2/dims.dh^2;
            % i-1,j
            AA(((j-1)*dims.nz+i-1)*3+1+1,((j-1)*dims.nz+i-2)*3+1+1)=C(3,3,i,j)/2/gz(i,j)^3/dims.dh*gz_z(i,j)+C(3,3,i,j)/gz(i,j)^2/dims.dh^2;
            % p
            % i+1,j
            AA(((j-1)*dims.nz+i-1)*3+1+1,((j-1)*dims.nz+i)*3+1+2)=-baralp3(i,j)/2/gz(i,j)/dims.dh;
            % i-1,j
            AA(((j-1)*dims.nz+i-1)*3+1+1,((j-1)*dims.nz+i-2)*3+1+2)=baralp3(i,j)/2/gz(i,j)/dims.dh;
            %% Eq. 3
            % ux
            % i,j+1
            AA(((j-1)*dims.nz+i-1)*3+1+2,((j)*dims.nz+i-1)*3+1)=-s3(i,j)*baralp1(i,j)/2/gx(i,j)/dims.dh;
            % i,j-1
            AA(((j-1)*dims.nz+i-1)*3+1+2,((j-2)*dims.nz+i-1)*3+1)=s3(i,j)*baralp1(i,j)/2/gx(i,j)/dims.dh;
            % uz
            % i+1,j
            AA(((j-1)*dims.nz+i-1)*3+1+2,((j-1)*dims.nz+i)*3+1+1)=-s1(i,j)*baralp3(i,j)/2/gz(i,j)/dims.dh;
            % i-1,j
            AA(((j-1)*dims.nz+i-1)*3+1+2,((j-1)*dims.nz+i-2)*3+1+1)=s1(i,j)*baralp3(i,j)/2/gz(i,j)/dims.dh;
            % p
            % i,j
            AA(((j-1)*dims.nz+i-1)*3+1+2,((j-1)*dims.nz+i-1)*3+1+2)=-bet(i,j)-2*s1(i,j)*the3(i,j)/gz(i,j)^2/dims.dh^2-2*s3(i,j)*the1(i,j)/gx(i,j)^2/dims.dh^2;
            % i,j+1
            AA(((j-1)*dims.nz+i-1)*3+1+2,((j)*dims.nz+i-1)*3+1+2)=s3(i,j)*the1(i,j)/2/gx(i,j)^3/dims.dh*gx_x(i,j)+s3(i,j)*the1(i,j)/gx(i,j)^2/dims.dh^2;
            % i,j-1
            AA(((j-1)*dims.nz+i-1)*3+1+2,((j-2)*dims.nz+i-1)*3+1+2)=-s3(i,j)*the1(i,j)/2/gx(i,j)^3/dims.dh*gx_x(i,j)+s3(i,j)*the1(i,j)/gx(i,j)^2/dims.dh^2;
            % i+1,j
            AA(((j-1)*dims.nz+i-1)*3+1+2,((j-1)*dims.nz+i)*3+1+2)=-s1(i,j)*the3(i,j)/2/gz(i,j)^3/dims.dh*gz_z(i,j)+s1(i,j)*the3(i,j)/gz(i,j)^2/dims.dh^2;
            % i-1,j
            AA(((j-1)*dims.nz+i-1)*3+1+2,((j-1)*dims.nz+i-2)*3+1+2)=s1(i,j)*the3(i,j)/2/gz(i,j)^3/dims.dh*gz_z(i,j)+s1(i,j)*the3(i,j)/gz(i,j)^2/dims.dh^2;
        end
    end
    %% upper left local
    for i=2:lp+1
        for j=2:lp+1
            %% Eq. 1
            % ux
            % i,j
            AA(((j-1)*dims.nz+i-1)*3+1,((j-1)*dims.nz+i-1)*3+1)=ome(l)^2*barrho1(i,j)-2*C(1,1,i,j)/gx(i,j)^2/dims.dh^2-2*C(4,4,i,j)/gz(i,j)^2/dims.dh^2;
            % i,j+1
            AA(((j-1)*dims.nz+i-1)*3+1,((j)*dims.nz+i-1)*3+1)=C(1,1,i,j)/2/gx(i,j)^3/dims.dh*gx_x(i,j)+C(1,1,i,j)/gx(i,j)^2/dims.dh^2;
            % i,j-1
            AA(((j-1)*dims.nz+i-1)*3+1,((j-2)*dims.nz+i-1)*3+1)=-C(1,1,i,j)/2/gx(i,j)^3/dims.dh*gx_x(i,j)+C(1,1,i,j)/gx(i,j)^2/dims.dh^2;
            % i+1,j
            AA(((j-1)*dims.nz+i-1)*3+1,((j-1)*dims.nz+i)*3+1)=C(4,4,i,j)/2/gz(i,j)^3/dims.dh*gz_z(i,j)+C(4,4,i,j)/gz(i,j)^2/dims.dh^2;
            % i-1,j
            AA(((j-1)*dims.nz+i-1)*3+1,((j-1)*dims.nz+i-2)*3+1)=-C(4,4,i,j)/2/gz(i,j)^3/dims.dh*gz_z(i,j)+C(4,4,i,j)/gz(i,j)^2/dims.dh^2;
            % uz
            % i+1,j+1
            AA(((j-1)*dims.nz+i-1)*3+1,((j)*dims.nz+i)*3+1+1)=(C(1,3,i,j)+C(4,4,i,j))/4/gx(i,j)/gz(i,j)/dims.dh^2;
            % i-1,j+1
            AA(((j-1)*dims.nz+i-1)*3+1,((j)*dims.nz+i-2)*3+1+1)=-(C(1,3,i,j)+C(4,4,i,j))/4/gx(i,j)/gz(i,j)/dims.dh^2;
            % i+1,j-1
            AA(((j-1)*dims.nz+i-1)*3+1,((j-2)*dims.nz+i)*3+1+1)=-(C(1,3,i,j)+C(4,4,i,j))/4/gx(i,j)/gz(i,j)/dims.dh^2;
            % i-1,j-1
            AA(((j-1)*dims.nz+i-1)*3+1,((j-2)*dims.nz+i-2)*3+1+1)=(C(1,3,i,j)+C(4,4,i,j))/4/gx(i,j)/gz(i,j)/dims.dh^2;
            % p
            % i,j+1
            AA(((j-1)*dims.nz+i-1)*3+1,((j)*dims.nz+i-1)*3+1+2)=-baralp1(i,j)/2/gx(i,j)/dims.dh;
            % i,j-1
            AA(((j-1)*dims.nz+i-1)*3+1,((j-2)*dims.nz+i-1)*3+1+2)=baralp1(i,j)/2/gx(i,j)/dims.dh;
            
            %% Eq. 2
            % ux
            % i+1,j+1
            AA(((j-1)*dims.nz+i-1)*3+1+1,((j)*dims.nz+i)*3+1)=(C(1,3,i,j)+C(4,4,i,j))/4/gx(i,j)/gz(i,j)/dims.dh^2;
            % i-1,j+1
            AA(((j-1)*dims.nz+i-1)*3+1+1,((j)*dims.nz+i-2)*3+1)=-(C(1,3,i,j)+C(4,4,i,j))/4/gx(i,j)/gz(i,j)/dims.dh^2;
            % i+1,j-1
            AA(((j-1)*dims.nz+i-1)*3+1+1,((j-2)*dims.nz+i)*3+1)=-(C(1,3,i,j)+C(4,4,i,j))/4/gx(i,j)/gz(i,j)/dims.dh^2;
            % i-1,j-1
            AA(((j-1)*dims.nz+i-1)*3+1+1,((j-2)*dims.nz+i-2)*3+1)=(C(1,3,i,j)+C(4,4,i,j))/4/gx(i,j)/gz(i,j)/dims.dh^2;
            % uz
            % i,j
            AA(((j-1)*dims.nz+i-1)*3+1+1,((j-1)*dims.nz+i-1)*3+1+1)=ome(l)^2*barrho3(i,j)-2*C(4,4,i,j)/gx(i,j)^2/dims.dh^2-2*C(3,3,i,j)/gz(i,j)^2/dims.dh^2;
            % i,j+1
            AA(((j-1)*dims.nz+i-1)*3+1+1,((j)*dims.nz+i-1)*3+1+1)=C(4,4,i,j)/2/gx(i,j)^3/dims.dh*gx_x(i,j)+C(4,4,i,j)/gx(i,j)^2/dims.dh^2;
            % i,j-1
            AA(((j-1)*dims.nz+i-1)*3+1+1,((j-2)*dims.nz+i-1)*3+1+1)=-C(4,4,i,j)/2/gx(i,j)^3/dims.dh*gx_x(i,j)+C(4,4,i,j)/gx(i,j)^2/dims.dh^2;
            % i+1,j
            AA(((j-1)*dims.nz+i-1)*3+1+1,((j-1)*dims.nz+i)*3+1+1)=C(3,3,i,j)/2/gz(i,j)^3/dims.dh*gz_z(i,j)+C(3,3,i,j)/gz(i,j)^2/dims.dh^2;
            % i-1,j
            AA(((j-1)*dims.nz+i-1)*3+1+1,((j-1)*dims.nz+i-2)*3+1+1)=-C(3,3,i,j)/2/gz(i,j)^3/dims.dh*gz_z(i,j)+C(3,3,i,j)/gz(i,j)^2/dims.dh^2;
            % p
            % i+1,j
            AA(((j-1)*dims.nz+i-1)*3+1+1,((j-1)*dims.nz+i)*3+1+2)=-baralp3(i,j)/2/gz(i,j)/dims.dh;
            % i-1,j
            AA(((j-1)*dims.nz+i-1)*3+1+1,((j-1)*dims.nz+i-2)*3+1+2)=baralp3(i,j)/2/gz(i,j)/dims.dh;
            %% Eq. 3
            % ux
            % i,j+1
            AA(((j-1)*dims.nz+i-1)*3+1+2,((j)*dims.nz+i-1)*3+1)=-s3(i,j)*baralp1(i,j)/2/gx(i,j)/dims.dh;
            % i,j-1
            AA(((j-1)*dims.nz+i-1)*3+1+2,((j-2)*dims.nz+i-1)*3+1)=s3(i,j)*baralp1(i,j)/2/gx(i,j)/dims.dh;
            % uz
            % i+1,j
            AA(((j-1)*dims.nz+i-1)*3+1+2,((j-1)*dims.nz+i)*3+1+1)=-s1(i,j)*baralp3(i,j)/2/gz(i,j)/dims.dh;
            % i-1,j
            AA(((j-1)*dims.nz+i-1)*3+1+2,((j-1)*dims.nz+i-2)*3+1+1)=s1(i,j)*baralp3(i,j)/2/gz(i,j)/dims.dh;
            % p
            % i,j
            AA(((j-1)*dims.nz+i-1)*3+1+2,((j-1)*dims.nz+i-1)*3+1+2)=-bet(i,j)-2*s1(i,j)*the3(i,j)/gz(i,j)^2/dims.dh^2-2*s3(i,j)*the1(i,j)/gx(i,j)^2/dims.dh^2;
            % i,j+1
            AA(((j-1)*dims.nz+i-1)*3+1+2,((j)*dims.nz+i-1)*3+1+2)=s3(i,j)*the1(i,j)/2/gx(i,j)^3/dims.dh*gx_x(i,j)+s3(i,j)*the1(i,j)/gx(i,j)^2/dims.dh^2;
            % i,j-1
            AA(((j-1)*dims.nz+i-1)*3+1+2,((j-2)*dims.nz+i-1)*3+1+2)=-s3(i,j)*the1(i,j)/2/gx(i,j)^3/dims.dh*gx_x(i,j)+s3(i,j)*the1(i,j)/gx(i,j)^2/dims.dh^2;
            % i+1,j
            AA(((j-1)*dims.nz+i-1)*3+1+2,((j-1)*dims.nz+i)*3+1+2)=s1(i,j)*the3(i,j)/2/gz(i,j)^3/dims.dh*gz_z(i,j)+s1(i,j)*the3(i,j)/gz(i,j)^2/dims.dh^2;
            % i-1,j
            AA(((j-1)*dims.nz+i-1)*3+1+2,((j-1)*dims.nz+i-2)*3+1+2)=-s1(i,j)*the3(i,j)/2/gz(i,j)^3/dims.dh*gz_z(i,j)+s1(i,j)*the3(i,j)/gz(i,j)^2/dims.dh^2;
        end
    end
    %% upper right local
    for i=2:lp+1
        for j=dims.nx-lp:dims.nx-1
            %% Eq. 1
            % ux
            % i,j
            AA(((j-1)*dims.nz+i-1)*3+1,((j-1)*dims.nz+i-1)*3+1)=ome(l)^2*barrho1(i,j)-2*C(1,1,i,j)/gx(i,j)^2/dims.dh^2-2*C(4,4,i,j)/gz(i,j)^2/dims.dh^2;
            % i,j+1
            AA(((j-1)*dims.nz+i-1)*3+1,((j)*dims.nz+i-1)*3+1)=-C(1,1,i,j)/2/gx(i,j)^3/dims.dh*gx_x(i,j)+C(1,1,i,j)/gx(i,j)^2/dims.dh^2;
            % i,j-1
            AA(((j-1)*dims.nz+i-1)*3+1,((j-2)*dims.nz+i-1)*3+1)=C(1,1,i,j)/2/gx(i,j)^3/dims.dh*gx_x(i,j)+C(1,1,i,j)/gx(i,j)^2/dims.dh^2;
            % i+1,j
            AA(((j-1)*dims.nz+i-1)*3+1,((j-1)*dims.nz+i)*3+1)=C(4,4,i,j)/2/gz(i,j)^3/dims.dh*gz_z(i,j)+C(4,4,i,j)/gz(i,j)^2/dims.dh^2;
            % i-1,j
            AA(((j-1)*dims.nz+i-1)*3+1,((j-1)*dims.nz+i-2)*3+1)=-C(4,4,i,j)/2/gz(i,j)^3/dims.dh*gz_z(i,j)+C(4,4,i,j)/gz(i,j)^2/dims.dh^2;
            % uz
            % i+1,j+1
            AA(((j-1)*dims.nz+i-1)*3+1,((j)*dims.nz+i)*3+1+1)=(C(1,3,i,j)+C(4,4,i,j))/4/gx(i,j)/gz(i,j)/dims.dh^2;
            % i-1,j+1
            AA(((j-1)*dims.nz+i-1)*3+1,((j)*dims.nz+i-2)*3+1+1)=-(C(1,3,i,j)+C(4,4,i,j))/4/gx(i,j)/gz(i,j)/dims.dh^2;
            % i+1,j-1
            AA(((j-1)*dims.nz+i-1)*3+1,((j-2)*dims.nz+i)*3+1+1)=-(C(1,3,i,j)+C(4,4,i,j))/4/gx(i,j)/gz(i,j)/dims.dh^2;
            % i-1,j-1
            AA(((j-1)*dims.nz+i-1)*3+1,((j-2)*dims.nz+i-2)*3+1+1)=(C(1,3,i,j)+C(4,4,i,j))/4/gx(i,j)/gz(i,j)/dims.dh^2;
            % p
            % i,j+1
            AA(((j-1)*dims.nz+i-1)*3+1,((j)*dims.nz+i-1)*3+1+2)=-baralp1(i,j)/2/gx(i,j)/dims.dh;
            % i,j-1
            AA(((j-1)*dims.nz+i-1)*3+1,((j-2)*dims.nz+i-1)*3+1+2)=baralp1(i,j)/2/gx(i,j)/dims.dh;
            
            %% Eq. 2
            % ux
            % i+1,j+1
            AA(((j-1)*dims.nz+i-1)*3+1+1,((j)*dims.nz+i)*3+1)=(C(1,3,i,j)+C(4,4,i,j))/4/gx(i,j)/gz(i,j)/dims.dh^2;
            % i-1,j+1
            AA(((j-1)*dims.nz+i-1)*3+1+1,((j)*dims.nz+i-2)*3+1)=-(C(1,3,i,j)+C(4,4,i,j))/4/gx(i,j)/gz(i,j)/dims.dh^2;
            % i+1,j-1
            AA(((j-1)*dims.nz+i-1)*3+1+1,((j-2)*dims.nz+i)*3+1)=-(C(1,3,i,j)+C(4,4,i,j))/4/gx(i,j)/gz(i,j)/dims.dh^2;
            % i-1,j-1
            AA(((j-1)*dims.nz+i-1)*3+1+1,((j-2)*dims.nz+i-2)*3+1)=(C(1,3,i,j)+C(4,4,i,j))/4/gx(i,j)/gz(i,j)/dims.dh^2;
            % uz
            % i,j
            AA(((j-1)*dims.nz+i-1)*3+1+1,((j-1)*dims.nz+i-1)*3+1+1)=ome(l)^2*barrho3(i,j)-2*C(4,4,i,j)/gx(i,j)^2/dims.dh^2-2*C(3,3,i,j)/gz(i,j)^2/dims.dh^2;
            % i,j+1
            AA(((j-1)*dims.nz+i-1)*3+1+1,((j)*dims.nz+i-1)*3+1+1)=-C(4,4,i,j)/2/gx(i,j)^3/dims.dh*gx_x(i,j)+C(4,4,i,j)/gx(i,j)^2/dims.dh^2;
            % i,j-1
            AA(((j-1)*dims.nz+i-1)*3+1+1,((j-2)*dims.nz+i-1)*3+1+1)=C(4,4,i,j)/2/gx(i,j)^3/dims.dh*gx_x(i,j)+C(4,4,i,j)/gx(i,j)^2/dims.dh^2;
            % i+1,j
            AA(((j-1)*dims.nz+i-1)*3+1+1,((j-1)*dims.nz+i)*3+1+1)=C(3,3,i,j)/2/gz(i,j)^3/dims.dh*gz_z(i,j)+C(3,3,i,j)/gz(i,j)^2/dims.dh^2;
            % i-1,j
            AA(((j-1)*dims.nz+i-1)*3+1+1,((j-1)*dims.nz+i-2)*3+1+1)=-C(3,3,i,j)/2/gz(i,j)^3/dims.dh*gz_z(i,j)+C(3,3,i,j)/gz(i,j)^2/dims.dh^2;
            % p
            % i+1,j
            AA(((j-1)*dims.nz+i-1)*3+1+1,((j-1)*dims.nz+i)*3+1+2)=-baralp3(i,j)/2/gz(i,j)/dims.dh;
            % i-1,j
            AA(((j-1)*dims.nz+i-1)*3+1+1,((j-1)*dims.nz+i-2)*3+1+2)=baralp3(i,j)/2/gz(i,j)/dims.dh;
            %% Eq. 3
            % ux
            % i,j+1
            AA(((j-1)*dims.nz+i-1)*3+1+2,((j)*dims.nz+i-1)*3+1)=-s3(i,j)*baralp1(i,j)/2/gx(i,j)/dims.dh;
            % i,j-1
            AA(((j-1)*dims.nz+i-1)*3+1+2,((j-2)*dims.nz+i-1)*3+1)=s3(i,j)*baralp1(i,j)/2/gx(i,j)/dims.dh;
            % uz
            % i+1,j
            AA(((j-1)*dims.nz+i-1)*3+1+2,((j-1)*dims.nz+i)*3+1+1)=-s1(i,j)*baralp3(i,j)/2/gz(i,j)/dims.dh;
            % i-1,j
            AA(((j-1)*dims.nz+i-1)*3+1+2,((j-1)*dims.nz+i-2)*3+1+1)=s1(i,j)*baralp3(i,j)/2/gz(i,j)/dims.dh;
            % p
            % i,j
            AA(((j-1)*dims.nz+i-1)*3+1+2,((j-1)*dims.nz+i-1)*3+1+2)=-bet(i,j)-2*s1(i,j)*the3(i,j)/gz(i,j)^2/dims.dh^2-2*s3(i,j)*the1(i,j)/gx(i,j)^2/dims.dh^2;
            % i,j+1
            AA(((j-1)*dims.nz+i-1)*3+1+2,((j)*dims.nz+i-1)*3+1+2)=-s3(i,j)*the1(i,j)/2/gx(i,j)^3/dims.dh*gx_x(i,j)+s3(i,j)*the1(i,j)/gx(i,j)^2/dims.dh^2;
            % i,j-1
            AA(((j-1)*dims.nz+i-1)*3+1+2,((j-2)*dims.nz+i-1)*3+1+2)=s3(i,j)*the1(i,j)/2/gx(i,j)^3/dims.dh*gx_x(i,j)+s3(i,j)*the1(i,j)/gx(i,j)^2/dims.dh^2;
            % i+1,j
            AA(((j-1)*dims.nz+i-1)*3+1+2,((j-1)*dims.nz+i)*3+1+2)=s1(i,j)*the3(i,j)/2/gz(i,j)^3/dims.dh*gz_z(i,j)+s1(i,j)*the3(i,j)/gz(i,j)^2/dims.dh^2;
            % i-1,j
            AA(((j-1)*dims.nz+i-1)*3+1+2,((j-1)*dims.nz+i-2)*3+1+2)=-s1(i,j)*the3(i,j)/2/gz(i,j)^3/dims.dh*gz_z(i,j)+s1(i,j)*the3(i,j)/gz(i,j)^2/dims.dh^2;
        end
    end
    %% lower left local
    for i=dims.nz-lp:dims.nz-1
        for j=2:lp+1
            %% Eq. 1
            % ux
            % i,j
            AA(((j-1)*dims.nz+i-1)*3+1,((j-1)*dims.nz+i-1)*3+1)=ome(l)^2*barrho1(i,j)-2*C(1,1,i,j)/gx(i,j)^2/dims.dh^2-2*C(4,4,i,j)/gz(i,j)^2/dims.dh^2;
            % i,j+1
            AA(((j-1)*dims.nz+i-1)*3+1,((j)*dims.nz+i-1)*3+1)=C(1,1,i,j)/2/gx(i,j)^3/dims.dh*gx_x(i,j)+C(1,1,i,j)/gx(i,j)^2/dims.dh^2;
            % i,j-1
            AA(((j-1)*dims.nz+i-1)*3+1,((j-2)*dims.nz+i-1)*3+1)=-C(1,1,i,j)/2/gx(i,j)^3/dims.dh*gx_x(i,j)+C(1,1,i,j)/gx(i,j)^2/dims.dh^2;
            % i+1,j
            AA(((j-1)*dims.nz+i-1)*3+1,((j-1)*dims.nz+i)*3+1)=-C(4,4,i,j)/2/gz(i,j)^3/dims.dh*gz_z(i,j)+C(4,4,i,j)/gz(i,j)^2/dims.dh^2;
            % i-1,j
            AA(((j-1)*dims.nz+i-1)*3+1,((j-1)*dims.nz+i-2)*3+1)=C(4,4,i,j)/2/gz(i,j)^3/dims.dh*gz_z(i,j)+C(4,4,i,j)/gz(i,j)^2/dims.dh^2;
            % uz
            % i+1,j+1
            AA(((j-1)*dims.nz+i-1)*3+1,((j)*dims.nz+i)*3+1+1)=(C(1,3,i,j)+C(4,4,i,j))/4/gx(i,j)/gz(i,j)/dims.dh^2;
            % i-1,j+1
            AA(((j-1)*dims.nz+i-1)*3+1,((j)*dims.nz+i-2)*3+1+1)=-(C(1,3,i,j)+C(4,4,i,j))/4/gx(i,j)/gz(i,j)/dims.dh^2;
            % i+1,j-1
            AA(((j-1)*dims.nz+i-1)*3+1,((j-2)*dims.nz+i)*3+1+1)=-(C(1,3,i,j)+C(4,4,i,j))/4/gx(i,j)/gz(i,j)/dims.dh^2;
            % i-1,j-1
            AA(((j-1)*dims.nz+i-1)*3+1,((j-2)*dims.nz+i-2)*3+1+1)=(C(1,3,i,j)+C(4,4,i,j))/4/gx(i,j)/gz(i,j)/dims.dh^2;
            % p
            % i,j+1
            AA(((j-1)*dims.nz+i-1)*3+1,((j)*dims.nz+i-1)*3+1+2)=-baralp1(i,j)/2/gx(i,j)/dims.dh;
            % i,j-1
            AA(((j-1)*dims.nz+i-1)*3+1,((j-2)*dims.nz+i-1)*3+1+2)=baralp1(i,j)/2/gx(i,j)/dims.dh;
            
            %% Eq. 2
            % ux
            % i+1,j+1
            AA(((j-1)*dims.nz+i-1)*3+1+1,((j)*dims.nz+i)*3+1)=(C(1,3,i,j)+C(4,4,i,j))/4/gx(i,j)/gz(i,j)/dims.dh^2;
            % i-1,j+1
            AA(((j-1)*dims.nz+i-1)*3+1+1,((j)*dims.nz+i-2)*3+1)=-(C(1,3,i,j)+C(4,4,i,j))/4/gx(i,j)/gz(i,j)/dims.dh^2;
            % i+1,j-1
            AA(((j-1)*dims.nz+i-1)*3+1+1,((j-2)*dims.nz+i)*3+1)=-(C(1,3,i,j)+C(4,4,i,j))/4/gx(i,j)/gz(i,j)/dims.dh^2;
            % i-1,j-1
            AA(((j-1)*dims.nz+i-1)*3+1+1,((j-2)*dims.nz+i-2)*3+1)=(C(1,3,i,j)+C(4,4,i,j))/4/gx(i,j)/gz(i,j)/dims.dh^2;
            % uz
            % i,j
            AA(((j-1)*dims.nz+i-1)*3+1+1,((j-1)*dims.nz+i-1)*3+1+1)=ome(l)^2*barrho3(i,j)-2*C(4,4,i,j)/gx(i,j)^2/dims.dh^2-2*C(3,3,i,j)/gz(i,j)^2/dims.dh^2;
            % i,j+1
            AA(((j-1)*dims.nz+i-1)*3+1+1,((j)*dims.nz+i-1)*3+1+1)=C(4,4,i,j)/2/gx(i,j)^3/dims.dh*gx_x(i,j)+C(4,4,i,j)/gx(i,j)^2/dims.dh^2;
            % i,j-1
            AA(((j-1)*dims.nz+i-1)*3+1+1,((j-2)*dims.nz+i-1)*3+1+1)=-C(4,4,i,j)/2/gx(i,j)^3/dims.dh*gx_x(i,j)+C(4,4,i,j)/gx(i,j)^2/dims.dh^2;
            % i+1,j
            AA(((j-1)*dims.nz+i-1)*3+1+1,((j-1)*dims.nz+i)*3+1+1)=-C(3,3,i,j)/2/gz(i,j)^3/dims.dh*gz_z(i,j)+C(3,3,i,j)/gz(i,j)^2/dims.dh^2;
            % i-1,j
            AA(((j-1)*dims.nz+i-1)*3+1+1,((j-1)*dims.nz+i-2)*3+1+1)=C(3,3,i,j)/2/gz(i,j)^3/dims.dh*gz_z(i,j)+C(3,3,i,j)/gz(i,j)^2/dims.dh^2;
            % p
            % i+1,j
            AA(((j-1)*dims.nz+i-1)*3+1+1,((j-1)*dims.nz+i)*3+1+2)=-baralp3(i,j)/2/gz(i,j)/dims.dh;
            % i-1,j
            AA(((j-1)*dims.nz+i-1)*3+1+1,((j-1)*dims.nz+i-2)*3+1+2)=baralp3(i,j)/2/gz(i,j)/dims.dh;
            %% Eq. 3
            % ux
            % i,j+1
            AA(((j-1)*dims.nz+i-1)*3+1+2,((j)*dims.nz+i-1)*3+1)=-s3(i,j)*baralp1(i,j)/2/gx(i,j)/dims.dh;
            % i,j-1
            AA(((j-1)*dims.nz+i-1)*3+1+2,((j-2)*dims.nz+i-1)*3+1)=s3(i,j)*baralp1(i,j)/2/gx(i,j)/dims.dh;
            % uz
            % i+1,j
            AA(((j-1)*dims.nz+i-1)*3+1+2,((j-1)*dims.nz+i)*3+1+1)=-s1(i,j)*baralp3(i,j)/2/gz(i,j)/dims.dh;
            % i-1,j
            AA(((j-1)*dims.nz+i-1)*3+1+2,((j-1)*dims.nz+i-2)*3+1+1)=s1(i,j)*baralp3(i,j)/2/gz(i,j)/dims.dh;
            % p
            % i,j
            AA(((j-1)*dims.nz+i-1)*3+1+2,((j-1)*dims.nz+i-1)*3+1+2)=-bet(i,j)-2*s1(i,j)*the3(i,j)/gz(i,j)^2/dims.dh^2-2*s3(i,j)*the1(i,j)/gx(i,j)^2/dims.dh^2;
            % i,j+1
            AA(((j-1)*dims.nz+i-1)*3+1+2,((j)*dims.nz+i-1)*3+1+2)=s3(i,j)*the1(i,j)/2/gx(i,j)^3/dims.dh*gx_x(i,j)+s3(i,j)*the1(i,j)/gx(i,j)^2/dims.dh^2;
            % i,j-1
            AA(((j-1)*dims.nz+i-1)*3+1+2,((j-2)*dims.nz+i-1)*3+1+2)=-s3(i,j)*the1(i,j)/2/gx(i,j)^3/dims.dh*gx_x(i,j)+s3(i,j)*the1(i,j)/gx(i,j)^2/dims.dh^2;
            % i+1,j
            AA(((j-1)*dims.nz+i-1)*3+1+2,((j-1)*dims.nz+i)*3+1+2)=-s1(i,j)*the3(i,j)/2/gz(i,j)^3/dims.dh*gz_z(i,j)+s1(i,j)*the3(i,j)/gz(i,j)^2/dims.dh^2;
            % i-1,j
            AA(((j-1)*dims.nz+i-1)*3+1+2,((j-1)*dims.nz+i-2)*3+1+2)=s1(i,j)*the3(i,j)/2/gz(i,j)^3/dims.dh*gz_z(i,j)+s1(i,j)*the3(i,j)/gz(i,j)^2/dims.dh^2;
        end
    end
    %% absorbing boundary ux
    % B1.2
    i=1;
    for j=2:dims.nx-1
        AA(3*((j-1)*dims.nz+i-1)+1,3*((j-1)*dims.nz+i-1)+1)=-1/dims.dh-1i*ome(l)/(vb);
        %AA((j-1)*dims.nz+i,(j)*dims.nz+i)=1/dims.dh^2;
        %AA((j-1)*dims.nz+i,(j-2)*dims.nz+i)=1/dims.dh^2;
        %AA((j-1)*dims.nz+i,(j-1)*dims.nz+i-1)=1/dims.dh^2;
        AA(3*((j-1)*dims.nz+i-1)+1,3*((j-1)*dims.nz+i)+1)=1/dims.dh;
    end
    % B1.3
    i=dims.nz;
    for j=2:dims.nx-1
        AA(3*((j-1)*dims.nz+i-1)+1,3*((j-1)*dims.nz+i-1)+1)=-1/dims.dh-1i*ome(l)/(vb);
        %AA((j-1)*dims.nz+i,(j)*dims.nz+i)=1/dims.dh;
        %AA((j-1)*dims.nz+i,(j-2)*dims.nz+i)=1/dims.dh;
        AA(3*((j-1)*dims.nz+i-1)+1,3*((j-1)*dims.nz+i-2)+1)=1/dims.dh;
        %AA((j-1)*dims.nz+i,(j-1)*dims.nz+i+1)=1/dims.dh;
    end
    % B1.4
    j=dims.nx;
    for i=2:dims.nz-1
        AA(3*((j-1)*dims.nz+i-1)+1,3*((j-1)*dims.nz+i-1)+1)=-1/dims.dh-1i*ome(l)/(vb);
        %AA((j-1)*dims.nz+i,(j)*dims.nz+i)=1/dims.dh;
        AA(3*((j-1)*dims.nz+i-1)+1,3*((j-2)*dims.nz+i-1)+1)=1/dims.dh;
        %AA((j-1)*dims.nz+i,(j-1)*dims.nz+i-1)=1/dims.dh;
        %AA((j-1)*dims.nz+i,(j-1)*dims.nz+i+1)=1/dims.dh;
    end
    % B1.5
    j=1;
    for i=2:dims.nz-1
        AA(3*((j-1)*dims.nz+i-1)+1,3*((j-1)*dims.nz+i-1)+1)=-1/dims.dh-1i*ome(l)/(vb);
        AA(3*((j-1)*dims.nz+i-1)+1,3*((j)*dims.nz+i-1)+1)=1/dims.dh;
        %AA((j-1)*dims.nz+i,(j-2)*dims.nz+i)=1/dims.dh;
        %AA((j-1)*dims.nz+i,(j-1)*dims.nz+i-1)=1/dims.dh;
        %AA((j-1)*dims.nz+i,(j-1)*dims.nz+i+1)=1/dims.dh;
    end
    % B1.6
    j=dims.nx;
    i=1;
    AA(3*((j-1)*dims.nz+i-1)+1,3*((j-1)*dims.nz+i-1)+1)=-2/dims.dh-2i*ome(l)/(vb);
    %AA((j-1)*dims.nz+i,(j)*dims.nz+i)=1/dims.dh;
    AA(3*((j-1)*dims.nz+i-1)+1,3*((j-1)*dims.nz+i)+1)=1/dims.dh;
    %AA((j-1)*dims.nz+i,(j-1)*dims.nz+i-1)=1/dims.dh;
    AA(3*((j-1)*dims.nz+i-1)+1,3*((j-2)*dims.nz+i-1)+1)=1/dims.dh;
    % B1.7
    j=1;
    i=1;
    AA(3*((j-1)*dims.nz+i-1)+1,3*((j-1)*dims.nz+i-1)+1)=-2/dims.dh-2i*ome(l)/(vb);
    AA(3*((j-1)*dims.nz+i-1)+1,3*((j-1)*dims.nz+i)+1)=1/dims.dh;
    %AA((j-1)*dims.nz+i,(j-2)*dims.nz+i)=1/dims.dh;
    %AA((j-1)*dims.nz+i,(j-1)*dims.nz+i-1)=1/dims.dh;
    AA(3*((j-1)*dims.nz+i-1)+1,3*((j)*dims.nz+i-1)+1)=1/dims.dh;
    % B1.8
    i=dims.nz;
    j=dims.nx;
    AA(3*((j-1)*dims.nz+i-1)+1,3*((j-1)*dims.nz+i-1)+1)=-2/dims.dh-2i*ome(l)/(vb);
    %AA((j-1)*dims.nz+i,(j)*dims.nz+i)=1/dims.dh;
    AA(3*((j-1)*dims.nz+i-1)+1,3*((j-1)*dims.nz+i-2)+1)=1/dims.dh;
    AA(3*((j-1)*dims.nz+i-1)+1,3*((j-2)*dims.nz+i-1)+1)=1/dims.dh;
    %AA((j-1)*dims.nz+i,(j-1)*dims.nz+i+1)=1/dims.dh;
    % B1.9
    j=1;
    i=dims.nz;
    AA(3*((j-1)*dims.nz+i-1)+1,3*((j-1)*dims.nz+i-1)+1)=-2/dims.dh-2i*ome(l)/(vb);
    AA(3*((j-1)*dims.nz+i-1)+1,3*((j-1)*dims.nz+i-2)+1)=1/dims.dh;
    %AA((j-1)*dims.nz+i,(j-2)*dims.nz+i)=1/dims.dh;
    AA(3*((j-1)*dims.nz+i-1)+1,3*((j)*dims.nz+i-1)+1)=1/dims.dh;
    %AA((j-1)*dims.nz+i,(j-1)*dims.nz+i+1)=1/dims.dh;
    %% absorbing boundary uz
    % B1.2
    i=1;
    for j=2:dims.nx-1
        AA(3*((j-1)*dims.nz+i-1)+1+1,3*((j-1)*dims.nz+i-1)+1+1)=-1/dims.dh-1i*ome(l)/(vb);
        %AA((j-1)*dims.nz+i,(j)*dims.nz+i)=1/dims.dh^2;
        %AA((j-1)*dims.nz+i,(j-2)*dims.nz+i)=1/dims.dh^2;
        %AA((j-1)*dims.nz+i,(j-1)*dims.nz+i-1)=1/dims.dh^2;
        AA(3*((j-1)*dims.nz+i-1)+1+1,3*((j-1)*dims.nz+i)+1+1)=1/dims.dh;
    end
    % B1.3
    i=dims.nz;
    for j=2:dims.nx-1
        AA(3*((j-1)*dims.nz+i-1)+1+1,3*((j-1)*dims.nz+i-1)+1+1)=-1/dims.dh-1i*ome(l)/(vb);
        %AA((j-1)*dims.nz+i,(j)*dims.nz+i)=1/dims.dh;
        %AA((j-1)*dims.nz+i,(j-2)*dims.nz+i)=1/dims.dh;
        AA(3*((j-1)*dims.nz+i-1)+1+1,3*((j-1)*dims.nz+i-2)+1+1)=1/dims.dh;
        %AA((j-1)*dims.nz+i,(j-1)*dims.nz+i+1)=1/dims.dh;
    end
    % B1.4
    j=dims.nx;
    for i=2:dims.nz-1
        AA(3*((j-1)*dims.nz+i-1)+1+1,3*((j-1)*dims.nz+i-1)+1+1)=-1/dims.dh-1i*ome(l)/(vb);
        %AA((j-1)*dims.nz+i,(j)*dims.nz+i)=1/dims.dh;
        AA(3*((j-1)*dims.nz+i-1)+1+1,3*((j-2)*dims.nz+i-1)+1+1)=1/dims.dh;
        %AA((j-1)*dims.nz+i,(j-1)*dims.nz+i-1)=1/dims.dh;
        %AA((j-1)*dims.nz+i,(j-1)*dims.nz+i+1)=1/dims.dh;
    end
    % B1.5
    j=1;
    for i=2:dims.nz-1
        AA(3*((j-1)*dims.nz+i-1)+1+1,3*((j-1)*dims.nz+i-1)+1+1)=-1/dims.dh-1i*ome(l)/(vb);
        AA(3*((j-1)*dims.nz+i-1)+1+1,3*((j)*dims.nz+i-1)+1+1)=1/dims.dh;
        %AA((j-1)*dims.nz+i,(j-2)*dims.nz+i)=1/dims.dh;
        %AA((j-1)*dims.nz+i,(j-1)*dims.nz+i-1)=1/dims.dh;
        %AA((j-1)*dims.nz+i,(j-1)*dims.nz+i+1)=1/dims.dh;
    end
    % B1.6
    j=dims.nx;
    i=1;
    AA(3*((j-1)*dims.nz+i-1)+1+1,3*((j-1)*dims.nz+i-1)+1+1)=-2/dims.dh-2i*ome(l)/(vb);
    %AA((j-1)*dims.nz+i,(j)*dims.nz+i)=1/dims.dh;
    AA(3*((j-1)*dims.nz+i-1)+1+1,3*((j-1)*dims.nz+i)+1+1)=1/dims.dh;
    %AA((j-1)*dims.nz+i,(j-1)*dims.nz+i-1)=1/dims.dh;
    AA(3*((j-1)*dims.nz+i-1)+1+1,3*((j-2)*dims.nz+i-1)+1+1)=1/dims.dh;
    % B1.7
    j=1;
    i=1;
    AA(3*((j-1)*dims.nz+i-1)+1+1,3*((j-1)*dims.nz+i-1)+1+1)=-2/dims.dh-2i*ome(l)/(vb);
    AA(3*((j-1)*dims.nz+i-1)+1+1,3*((j-1)*dims.nz+i)+1+1)=1/dims.dh;
    %AA((j-1)*dims.nz+i,(j-2)*dims.nz+i)=1/dims.dh;
    %AA((j-1)*dims.nz+i,(j-1)*dims.nz+i-1)=1/dims.dh;
    AA(3*((j-1)*dims.nz+i-1)+1+1,3*((j)*dims.nz+i-1)+1+1)=1/dims.dh;
    % B1.8
    i=dims.nz;
    j=dims.nx;
    AA(3*((j-1)*dims.nz+i-1)+1+1,3*((j-1)*dims.nz+i-1)+1+1)=-2/dims.dh-2i*ome(l)/(vb);
    %AA((j-1)*dims.nz+i,(j)*dims.nz+i)=1/dims.dh;
    AA(3*((j-1)*dims.nz+i-1)+1+1,3*((j-1)*dims.nz+i-2)+1+1)=1/dims.dh;
    AA(3*((j-1)*dims.nz+i-1)+1+1,3*((j-2)*dims.nz+i-1)+1+1)=1/dims.dh;
    %AA((j-1)*dims.nz+i,(j-1)*dims.nz+i+1)=1/dims.dh;
    % B1.9
    j=1;
    i=dims.nz;
    AA(3*((j-1)*dims.nz+i-1)+1+1,3*((j-1)*dims.nz+i-1)+1+1)=-2/dims.dh-2i*ome(l)/(vb);
    AA(3*((j-1)*dims.nz+i-1)+1+1,3*((j-1)*dims.nz+i-2)+1+1)=1/dims.dh;
    %AA((j-1)*dims.nz+i,(j-2)*dims.nz+i)=1/dims.dh;
    AA(3*((j-1)*dims.nz+i-1)+1+1,3*((j)*dims.nz+i-1)+1+1)=1/dims.dh;
    %AA((j-1)*dims.nz+i,(j-1)*dims.nz+i+1)=1/dims.dh;
    %% p
    % B1.2
    i=1;
    for j=2:dims.nx-1
        AA(3*((j-1)*dims.nz+i-1)+1+2,3*((j-1)*dims.nz+i-1)+1+2)=-1/dims.dh-1i*ome(l)/(vb);
        %AA((j-1)*dims.nz+i,(j)*dims.nz+i)=1/dims.dh^2;
        %AA((j-1)*dims.nz+i,(j-2)*dims.nz+i)=1/dims.dh^2;
        %AA((j-1)*dims.nz+i,(j-1)*dims.nz+i-1)=1/dims.dh^2;
        AA(3*((j-1)*dims.nz+i-1)+1+2,3*((j-1)*dims.nz+i)+1+2)=1/dims.dh;
    end
    % B1.3
    i=dims.nz;
    for j=2:dims.nx-1
        AA(3*((j-1)*dims.nz+i-1)+1+2,3*((j-1)*dims.nz+i-1)+1+2)=-1/dims.dh-1i*ome(l)/(vb);
        %AA((j-1)*dims.nz+i,(j)*dims.nz+i)=1/dims.dh;
        %AA((j-1)*dims.nz+i,(j-2)*dims.nz+i)=1/dims.dh;
        AA(3*((j-1)*dims.nz+i-1)+1+2,3*((j-1)*dims.nz+i-2)+1+2)=1/dims.dh;
        %AA((j-1)*dims.nz+i,(j-1)*dims.nz+i+1)=1/dims.dh;
    end
    % B1.4
    j=dims.nx;
    for i=2:dims.nz-1
        AA(3*((j-1)*dims.nz+i-1)+1+2,3*((j-1)*dims.nz+i-1)+1+2)=-1/dims.dh-1i*ome(l)/(vb);
        %AA((j-1)*dims.nz+i,(j)*dims.nz+i)=1/dims.dh;
        AA(3*((j-1)*dims.nz+i-1)+1+2,3*((j-2)*dims.nz+i-1)+1+2)=1/dims.dh;
        %AA((j-1)*dims.nz+i,(j-1)*dims.nz+i-1)=1/dims.dh;
        %AA((j-1)*dims.nz+i,(j-1)*dims.nz+i+1)=1/dims.dh;
    end
    % B1.5
    j=1;
    for i=2:dims.nz-1
        AA(3*((j-1)*dims.nz+i-1)+1+2,3*((j-1)*dims.nz+i-1)+1+2)=-1/dims.dh-1i*ome(l)/(vb);
        AA(3*((j-1)*dims.nz+i-1)+1+2,3*((j)*dims.nz+i-1)+1+2)=1/dims.dh;
        %AA((j-1)*dims.nz+i,(j-2)*dims.nz+i)=1/dims.dh;
        %AA((j-1)*dims.nz+i,(j-1)*dims.nz+i-1)=1/dims.dh;
        %AA((j-1)*dims.nz+i,(j-1)*dims.nz+i+1)=1/dims.dh;
    end
    % B1.6
    j=dims.nx;
    i=1;
    AA(3*((j-1)*dims.nz+i-1)+1+2,3*((j-1)*dims.nz+i-1)+1+2)=-2/dims.dh-2i*ome(l)/(vb);
    %AA((j-1)*dims.nz+i,(j)*dims.nz+i)=1/dims.dh;
    AA(3*((j-1)*dims.nz+i-1)+1+2,3*((j-1)*dims.nz+i)+1+2)=1/dims.dh;
    %AA((j-1)*dims.nz+i,(j-1)*dims.nz+i-1)=1/dims.dh;
    AA(3*((j-1)*dims.nz+i-1)+1+2,3*((j-2)*dims.nz+i-1)+1+2)=1/dims.dh;
    % B1.7
    j=1;
    i=1;
    AA(3*((j-1)*dims.nz+i-1)+1+2,3*((j-1)*dims.nz+i-1)+1+2)=-2/dims.dh-2i*ome(l)/(vb);
    AA(3*((j-1)*dims.nz+i-1)+1+2,3*((j-1)*dims.nz+i)+1+2)=1/dims.dh;
    %AA((j-1)*dims.nz+i,(j-2)*dims.nz+i)=1/dims.dh;
    %AA((j-1)*dims.nz+i,(j-1)*dims.nz+i-1)=1/dims.dh;
    AA(3*((j-1)*dims.nz+i-1)+1+2,3*((j)*dims.nz+i-1)+1+2)=1/dims.dh;
    % B1.8
    i=dims.nz;
    j=dims.nx;
    AA(3*((j-1)*dims.nz+i-1)+1+2,3*((j-1)*dims.nz+i-1)+1+2)=-2/dims.dh-2i*ome(l)/(vb);
    %AA((j-1)*dims.nz+i,(j)*dims.nz+i)=1/dims.dh;
    AA(3*((j-1)*dims.nz+i-1)+1+2,3*((j-1)*dims.nz+i-2)+1+2)=1/dims.dh;
    AA(3*((j-1)*dims.nz+i-1)+1+2,3*((j-2)*dims.nz+i-1)+1+2)=1/dims.dh;
    %AA((j-1)*dims.nz+i,(j-1)*dims.nz+i+1)=1/dims.dh;
    % B1.9
    j=1;
    i=dims.nz;
    AA(3*((j-1)*dims.nz+i-1)+1+2,3*((j-1)*dims.nz+i-1)+1+2)=-2/dims.dh-2i*ome(l)/(vb);
    AA(3*((j-1)*dims.nz+i-1)+1+2,3*((j-1)*dims.nz+i-2)+1+2)=1/dims.dh;
    %AA((j-1)*dims.nz+i,(j-2)*dims.nz+i)=1/dims.dh;
    AA(3*((j-1)*dims.nz+i-1)+1+2,3*((j)*dims.nz+i-1)+1+2)=1/dims.dh;
    %AA((j-1)*dims.nz+i,(j-1)*dims.nz+i+1)=1/dims.dh;
    %% source
    s(3*((dims.sx-1)*dims.nz+dims.sz-1)+1)=dims.s1(l);
    s(3*((dims.sx-1)*dims.nz+dims.sz-1)+1+1)=dims.s2(l);
    s(3*((dims.sx-1)*dims.nz+dims.sz-1)+1+2)=dims.s3(l);
    u=AA\(-s);
    %% assemble ux,uz,p
    for i=1:length(u)
        if mod(i,3)==1
            ux(fix(i/3)+1+(l-1)*dims.nz*dims.nx)=u(i);
        end
        if mod(i,3)==2
            uz(fix(i/3)+1+(l-1)*dims.nz*dims.nx)=u(i);
        end
        if mod(i,3)==0
            p(fix(i/3)+(l-1)*dims.nz*dims.nx)=u(i);
        end
    end
    fprintf('l=%f/%f\n',l,length(ome));
    toc;
end
%%
uxt=ifft(ux,dims.nt,3)*dims.nt;
uzt=ifft(uz,dims.nt,3)*dims.nt;
pt=ifft(p,dims.nt,3)*dims.nt;
recx=reshape(uxt(dims.rz(1),dims.rx,:),[length(dims.rx),dims.nt]).';
recz=reshape(uzt(dims.rz(1),dims.rx,:),[length(dims.rx),dims.nt]).';
recp=reshape(pt(dims.rz(1),dims.rx,:),[length(dims.rx),dims.nt]).';
%%
end
