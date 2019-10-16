%% Dimensions
close all;
load('kl.mat');
load('C_b.mat');

dims.dh=8; % Spacial grid step
dims.dt=10^-2; % [s]

dims.nt=250; % Amount of time steps
dims.ns=250;

b=zeros(4,1);
b(1)=40;
b(2)=b(1)+fix((kl.i(2)-kl.i(1))/dims.dh);
b(3)=fix((kl.i(3)-kl.i(2))/dims.dh)+b(2);
b(4)=fix((kl.i(4)-kl.i(3))/dims.dh)+b(3);
dims.nz=b(4)+40+10;
dims.nx=b(4)*2+40;
%% Model dimensions
dims.mz=b(1)+1:dims.nz-39;
dims.mx=40:dims.nx-39;
%% Source and source signals
dims.sz=b(1)+1;
dims.sx=fix(dims.nx/2);
sn=length(dims.sx);
singles=rickerWave(10,dims);
figure(1)
plot(singles);
%% Receiver locations
dims.rx=min(dims.mx):max(dims.mx);
dims.rz=min(dims.mz)*ones(1,length(dims.rx));
%% Assign physical parameters to poroelastic layer
vp03=kl.vp0(3);
vs03=kl.vs0(3);
rho3=kl.rho(3);
C33_3=vp03^2*rho3;
C44_3=vs03^2*rho3;
C_i=C_b; % C of modified model with an isotropic layer
C_i(3,3,3)=C33_3;
C_i(4,4,3)=C44_3;
C_i(1,1,3)=C_i(3,3,3);
C_i(2,2,3)=C_i(3,3,3);
C_i(5,5,3)=C_i(4,4,3);
C_i(6,6,3)=C_i(4,4,3);
C_i(1,2,3)=C_i(3,3,3)-2*C_i(4,4,3);
C_i(1,3,3)=C_i(1,2,3);
C_i(2,3,3)=C_i(1,2,3);



Ks=zeros(dims.nz,dims.nx);
Ks(:)=38*10^9;

rhos=zeros(dims.nz,dims.nx);
rhos(:)=2650;

Kf=zeros(dims.nz,dims.nx);
Kf(:)=8.64*10^8;

rhof=zeros(dims.nz,dims.nx);
rhof(:)=600;

eta=zeros(dims.nz,dims.nx);
eta(:)=10^-1;

phi=zeros(dims.nz,dims.nx);
phi(b(3)+1:b(4),:)=0.001;

k3=zeros(dims.nz,dims.nx);
k3(:)=10^-15;

k1=zeros(dims.nz,dims.nx);
k1(:)=10^-15;

C=zeros(6,6,dims.nz,dims.nx);

% block 1
C(1,1,1:b(2),:)=C_i(1,1,1);
C(1,2,1:b(2),:)=C_i(1,2,1);
C(1,3,1:b(2),:)=C_i(1,3,1);
C(3,3,1:b(2),:)=C_i(3,3,1);
C(4,4,1:b(2),:)=C_i(4,4,1);

% block 2
C(1,1,b(2)+1:b(3),:)=C_i(1,1,2);
C(1,2,b(2)+1:b(3),:)=C_i(1,2,2);
C(1,3,b(2)+1:b(3),:)=C_i(1,3,2);
C(3,3,b(2)+1:b(3),:)=C_i(3,3,2);
C(4,4,b(2)+1:b(3),:)=C_i(4,4,2);

% block 3
C(1,1,b(3)+1:b(4),:)=C_i(1,1,3);
C(1,2,b(3)+1:b(4),:)=C_i(1,2,3);
C(1,3,b(3)+1:b(4),:)=C_i(1,3,3);
C(3,3,b(3)+1:b(4),:)=C_i(3,3,3);
C(4,4,b(3)+1:b(4),:)=C_i(4,4,3);

% block 4
C(1,1,b(4)+1:end,:)=C_i(1,1,4);
C(1,2,b(4)+1:end,:)=C_i(1,2,4);
C(1,3,b(4)+1:end,:)=C_i(1,3,4);
C(3,3,b(4)+1:end,:)=C_i(3,3,4);
C(4,4,b(4)+1:end,:)=C_i(4,4,4);

R1=zeros(dims.nz,dims.nx);
R1(:)=10^-3;
R3=R1;

rhoa=zeros(dims.nz,dims.nx);
rhoa(:)=550;
%% source
fs=1/dims.dt;
L=dims.nt;
n=L;
f=fs*(0:(n/2))/n;
s=zeros([dims.nz*dims.nx,1]);

source_freq=fft(singles,n)/(n/2);
source_freq2=source_freq(1:n/2+1);

plot(f,abs(source_freq2));
xlabel('f [Hz]');
ylabel('|A| [m]');
%% check source
source_time=ifft(source_freq2,n,1)*n;
figure(2)
subplot(2,1,1)
plot(dims.dt:dims.dt:dims.dt*L,real(source_time));
subplot(2,1,2)
plot(dims.dt:dims.dt:dims.dt*L,real(singles));
hold off;
%%
s_diff=diff(abs(source_freq2));
s_diff=[s_diff(1);s_diff];
s_lim=find(abs(source_freq2)<.01*max(abs(source_freq2)) & s_diff<0);

s_lim2=s_lim(1);
f_range=1:s_lim2;
f2=f(f_range);
ome=2*pi*f2;
dims.s1=0*source_freq2(f_range);
dims.s2=1*source_freq2(f_range);
dims.s3=0*source_freq2(f_range);
%% boundary condition input
lp=40;
d0=-3*3500*log(10^-6)/2/40/dims.dh;
vb=3500;
%%
[ux,uz,p,uxt,uzt,pt,recx,recz,recp]=poroelastic_solver_PLM(dims,C,Ks,rhos,Kf,rhof,eta,phi,k1,k3,R1,R3,rhoa,ome,lp,d0,vb);
%% model
%{
C2=reshape(C(3,3,:,:),[size(C(3,3,:,:),3),size(C(3,3,:,:),4)]);
figure(7)
set(gcf,'position',[80,80,1000,700]);
imagesc(C2);
colorbar;
hold on;
ax=plot([min(dims.rx),max(dims.rx)],[dims.rz(1),dims.rz(1)],'color','red','linewidth',2);
hold on;
ax2=scatter(dims.sx,dims.sz,30,[0,1,0],'filled');
hold on;
%{
ax3=plot([lp+1,lp+1],[lp+1,dims.nz-lp-1],'color',[1,.3,1]);
hold on;
ax3=plot([dims.nx-lp-1,dims.nx-lp-1],[lp+1,dims.nz-lp-1],'color',[1,.3,1]);
hold on;
ax3=plot([lp+1,dims.nx-lp-1],[lp+1,lp+1],'color',[1,.3,1]);
hold on;
ax3=plot([lp+1,dims.nx-lp-1],[dims.nz-lp-1,dims.nz-lp-1],'color',[1,.3,1]);
axis on;
%}
xlabel({['x*' num2str(dims.dh) '[m]']});
ylabel({['z*' num2str(dims.dh) '[m]']});
title({['C_{33}']});
legend([ax2,ax],'source','receiver','Location',[0.48,0.03,0.005,0.005],'orientation','horizontal');
print(gcf,['C:\Users\Yi\OneDrive\master thesis\VTI_poroelasticity\sch_model'],'-dpng','-r600');
%{
aa=reshape(uzt(:,:,200),[size(uzt,1),size(uzt,2)]);
figure(8)
fu=imfuse(real(aa),C33,'falsecolor','colorchannels',[0,1,2]);
imagesc(fu);
%}
%}
%% frequency domain
for i=1:2:size(ux,3)
    figure(4)
    set(gcf,'Visible','on');
    set(gcf,'position',[80,80,1000,600]);
    
    subplot(2,2,1)
    imshow(real(ux(:,:,i)),1*[min(min(real(ux(:,:,i)))),max(max(real(ux(:,:,i))))]);
    colorbar;
    hold on;
    ax=plot([min(dims.rx),max(dims.rx)],[dims.rz(1),dims.rz(1)],'color','red','linewidth',2);
    hold on;
    ax2=scatter(dims.sx,dims.sz,30,[0,1,0],'filled');
    hold on;
    ax3=plot([lp+1,lp+1],[lp+1,dims.nz-lp-1],'color','blue');
    hold on;
    ax3=plot([dims.nx-lp-1,dims.nx-lp-1],[lp+1,dims.nz-lp-1],'color','blue');
    hold on;
    ax3=plot([lp+1,dims.nx-lp-1],[lp+1,lp+1],'color','blue');
    hold on;
    ax3=plot([lp+1,dims.nx-lp-1],[dims.nz-lp-1,dims.nz-lp-1],'color','blue');
    hold on;
    ax4=plot([min(dims.mx),max(dims.mx)],[b(2),b(2)],'color','yellow');
    hold on;
    ax4=plot([min(dims.mx),max(dims.mx)],[b(3),b(3)],'color','yellow');
    hold on;
    ax4=plot([min(dims.mx),max(dims.mx)],[b(4),b(4)],'color','yellow');
    axis on;
    xlabel({['x*' num2str(dims.dh) '[m]']});
    ylabel({['z*' num2str(dims.dh) '[m]']});
    title(sprintf('ome=%frad/s \nux',ome(i)));
    
    subplot(2,2,2)
    imshow(real(uz(:,:,i)),1*[min(min(real(uz(:,:,i)))),max(max(real(uz(:,:,i))))]);
    colorbar;
    hold on;
    ax=plot([min(dims.rx),max(dims.rx)],[dims.rz(1),dims.rz(1)],'color','red','linewidth',2);
    hold on;
    ax2=scatter(dims.sx,dims.sz,30,[0,1,0],'filled');
    hold on;
    ax3=plot([lp+1,lp+1],[lp+1,dims.nz-lp-1],'color','blue');
    hold on;
    ax3=plot([dims.nx-lp-1,dims.nx-lp-1],[lp+1,dims.nz-lp-1],'color','blue');
    hold on;
    ax3=plot([lp+1,dims.nx-lp-1],[lp+1,lp+1],'color','blue');
    hold on;
    ax3=plot([lp+1,dims.nx-lp-1],[dims.nz-lp-1,dims.nz-lp-1],'color','blue');
    hold on;
    ax4=plot([min(dims.mx),max(dims.mx)],[b(2),b(2)],'color','yellow');
    hold on;
    ax4=plot([min(dims.mx),max(dims.mx)],[b(3),b(3)],'color','yellow');
    hold on;
    ax4=plot([min(dims.mx),max(dims.mx)],[b(4),b(4)],'color','yellow');
    axis on;
    title('uz');
    xlabel({['x*' num2str(dims.dh) '[m]']});
    ylabel({['z*' num2str(dims.dh) '[m]']});
    
    subplot(2,2,3)
    imshow(real(p(:,:,i)),1*[min(min(real(p(:,:,i)))),max(max(real(p(:,:,i))))]);
    colorbar;
    hold on;
    ax=plot([min(dims.rx),max(dims.rx)],[dims.rz(1),dims.rz(1)],'color','red','linewidth',2);
    hold on;
    ax2=scatter(dims.sx,dims.sz,30,[0,1,0],'filled');
    hold on;
    ax3=plot([lp+1,lp+1],[lp+1,dims.nz-lp-1],'color','blue');
    hold on;
    ax3=plot([dims.nx-lp-1,dims.nx-lp-1],[lp+1,dims.nz-lp-1],'color','blue');
    hold on;
    ax3=plot([lp+1,dims.nx-lp-1],[lp+1,lp+1],'color','blue');
    hold on;
    ax3=plot([lp+1,dims.nx-lp-1],[dims.nz-lp-1,dims.nz-lp-1],'color','blue');
    hold on;
    ax4=plot([min(dims.mx),max(dims.mx)],[b(2),b(2)],'color','yellow');
    hold on;
    ax4=plot([min(dims.mx),max(dims.mx)],[b(3),b(3)],'color','yellow');
    hold on;
    ax4=plot([min(dims.mx),max(dims.mx)],[b(4),b(4)],'color','yellow');
    axis on;
    title('p');
    xlabel({['x*' num2str(dims.dh) '[m]']});
    ylabel({['z*' num2str(dims.dh) '[m]']});
    
    legend([ax2,ax,ax3,ax4],'source','receiver','PML boundary','block boundary','Location',[0.5,0.03,0.005,0.005],'orientation','horizontal');
    
    shg;
    %print(gcf,['C:\Users\Yi\OneDrive\master thesis\VTI_poroelasticity\oil_10hz\freq' num2str(i) '.png'],'-dpng','-r600');
end
hold off;
%% time domain
for i=1:5:120
    figure(4)
    set(gcf,'Visible','on');
    set(gcf,'position',[80,80,1000,600]);
    
    subplot(2,2,1)
    imshow(real(uxt(:,:,i)),.1*[min(min(real(uxt(:)))),max(max(real(uxt(:))))]);
    colorbar;
    %grid on;
    hold on;
    ax=plot([min(dims.rx),max(dims.rx)],[dims.rz(1),dims.rz(1)],'color','red','linewidth',2);
    hold on;
    ax2=scatter(dims.sx,dims.sz,30,[0,1,0],'filled');
    hold on;
    ax3=plot([lp+1,lp+1],[lp+1,dims.nz-lp-1],'color','blue');
    hold on;
    ax3=plot([dims.nx-lp-1,dims.nx-lp-1],[lp+1,dims.nz-lp-1],'color','blue');
    hold on;
    ax3=plot([lp+1,dims.nx-lp-1],[lp+1,lp+1],'color','blue');
    hold on;
    ax3=plot([lp+1,dims.nx-lp-1],[dims.nz-lp-1,dims.nz-lp-1],'color','blue');
    hold on;
    ax4=plot([min(dims.mx),max(dims.mx)],[b(2),b(2)],'color','yellow');
    hold on;
    ax4=plot([min(dims.mx),max(dims.mx)],[b(3),b(3)],'color','yellow');
    hold on;
    ax4=plot([min(dims.mx),max(dims.mx)],[b(4),b(4)],'color','yellow');
    axis on;
    xlabel({['x*' num2str(dims.dh) '[m]']});
    ylabel({['z*' num2str(dims.dh) '[m]']});
    title({['t=' num2str(i*dims.dt) 's'],['uxt']});
    hold on;
    
    subplot(2,2,2)
    imshow(real(uzt(:,:,i)),.1*[min(min(real(uzt(:)))),max(max(real(uzt(:))))]);
    colorbar;
    %grid on;
    hold on;
    ax=plot([min(dims.rx),max(dims.rx)],[dims.rz(1),dims.rz(1)],'color','red','linewidth',2);
    hold on;
    ax2=scatter(dims.sx,dims.sz,30,[0,1,0],'filled');
    hold on;
    ax3=plot([lp+1,lp+1],[lp+1,dims.nz-lp-1],'color','blue');
    hold on;
    ax3=plot([dims.nx-lp-1,dims.nx-lp-1],[lp+1,dims.nz-lp-1],'color','blue');
    hold on;
    ax3=plot([lp+1,dims.nx-lp-1],[lp+1,lp+1],'color','blue');
    hold on;
    ax3=plot([lp+1,dims.nx-lp-1],[dims.nz-lp-1,dims.nz-lp-1],'color','blue');
    hold on;
    ax4=plot([min(dims.mx),max(dims.mx)],[b(2),b(2)],'color','yellow');
    hold on;
    ax4=plot([min(dims.mx),max(dims.mx)],[b(3),b(3)],'color','yellow');
    hold on;
    ax4=plot([min(dims.mx),max(dims.mx)],[b(4),b(4)],'color','yellow');
    axis on;
    xlabel({['x*' num2str(dims.dh) '[m]']});
    ylabel({['z*' num2str(dims.dh) '[m]']});
    title({['uzt']});
    hold on;
    
    subplot(2,2,3)
    imshow(real(pt(:,:,i)),.1*[min(real(pt(:))),max(real(pt(:)))]);
    colorbar;
    %grid on;
    hold on;
    ax=plot([min(dims.rx),max(dims.rx)],[dims.rz(1),dims.rz(1)],'color','red','linewidth',2);
    hold on;
    ax2=scatter(dims.sx,dims.sz,30,[0,1,0],'filled');
    hold on;
    ax3=plot([lp+1,lp+1],[lp+1,dims.nz-lp-1],'color','blue');
    hold on;
    ax3=plot([dims.nx-lp-1,dims.nx-lp-1],[lp+1,dims.nz-lp-1],'color','blue');
    hold on;
    ax3=plot([lp+1,dims.nx-lp-1],[lp+1,lp+1],'color','blue');
    hold on;
    ax3=plot([lp+1,dims.nx-lp-1],[dims.nz-lp-1,dims.nz-lp-1],'color','blue');
    hold on;
    ax4=plot([min(dims.mx),max(dims.mx)],[b(2),b(2)],'color','yellow');
    hold on;
    ax4=plot([min(dims.mx),max(dims.mx)],[b(3),b(3)],'color','yellow');
    hold on;
    ax4=plot([min(dims.mx),max(dims.mx)],[b(4),b(4)],'color','yellow');
    axis on;
    xlabel({['x*' num2str(dims.dh) '[m]']});
    ylabel({['z*' num2str(dims.dh) '[m]']});
    title({['pt']});
    hold on;
    
    legend([ax2,ax,ax3,ax4],'source','receiver','PML boundary','block boundary','Location',[0.5,0.03,0.005,0.005],'orientation','horizontal');
    shg;
    %print(gcf,['C:\Users\Yi\OneDrive\master thesis\VTI_poroelasticity\gas_10hz\time' num2str(i) '.png'],'-dpng','-r600');
end
hold off;
%% recordings
for i=1:10:size(uzt,3)
    figure(5)
    set(gcf,'Visible','on');
    set(gcf,'position',[80,80,1500,600]);
    
    subplot(1,2,1)
    imagesc(real(uzt(:,:,i)),.1*[min(min(real(uzt(:)))),max(max(real(uzt(:))))]);
    colorbar;
    %grid on;
    hold on;
    ax=plot([min(dims.rx),max(dims.rx)],[dims.rz(1),dims.rz(1)],'color','red','linewidth',2);
    hold on;
    ax2=scatter(dims.sx,dims.sz,30,[0,1,0],'filled');
    hold on;
    ax3=plot([lp+1,lp+1],[lp+1,dims.nz-lp-1],'color','blue');
    hold on;
    ax3=plot([dims.nx-lp-1,dims.nx-lp-1],[lp+1,dims.nz-lp-1],'color','blue');
    hold on;
    ax3=plot([lp+1,dims.nx-lp-1],[lp+1,lp+1],'color','blue');
    hold on;
    ax3=plot([lp+1,dims.nx-lp-1],[dims.nz-lp-1,dims.nz-lp-1],'color','blue');
    hold on;
    ax4=plot([min(dims.mx),max(dims.mx)],[b(2),b(2)],'color','yellow');
    hold on;
    ax4=plot([min(dims.mx),max(dims.mx)],[b(3),b(3)],'color','yellow');
    hold on;
    ax4=plot([min(dims.mx),max(dims.mx)],[b(4),b(4)],'color','yellow');
    axis on;
    xlabel({['x*' num2str(dims.dh) '[m]']});
    ylabel({['z*' num2str(dims.dh) '[m]']});
    title({['t=' num2str(i*dims.dt) 's'],['uzt']});
    
    subplot(1,2,2)
    imagesc(real(recz(1:i,:)),.1*[min(real(recz(:))),max(real(recz(:)))]);
    colorbar;
    axis on;
    %grid on;
    ylim([0,250]);
    xlabel({['x*' num2str(dims.dh) '[m]']});
    ylabel({['t*' num2str(dims.dt) '[s]']});
    title({['recz']});
    
    legend([ax2,ax,ax3,ax4],'source','receiver','PML boundary','block boundary','Location',[0.5,0.03,0.005,0.005],'orientation','horizontal');
    hold on;
    shg;
    %print(gcf,['C:\Users\Yi\OneDrive\master thesis\VTI_poroelasticity\oil_10hz\recz\' num2str(i) '.png'],'-dpng','-r600');
end
%% oil
%{
save('oil_ux_phi_.001_10hz.mat','ux');
save('oil_uz_phi_.001_10hz.mat','uz');
save('oil_p_phi_.001_10hz.mat','p');
save('oil_uxt_phi_.001_10hz.mat','uxt');
save('oil_uzt_phi_.001_10hz.mat','uzt');
save('oil_pt_phi_.001_10hz.mat','pt');
save('oil_recx_phi_.001_10hz.mat','recx');
save('oil_recz_phi_.001_10hz.mat','recz');
%% water
Kf(:)=2.25*10^9;
rhof(:)=1000;
eta(:)=10^-3;
[ux,uz,p,uxt,uzt,pt,recx,recz,recp]= ... 
    poroelastic_solver_PLM(dims,C,Ks,rhos,Kf,rhof,eta,phi,k1,k3,R1,R3,rhoa,ome,lp,d0,vb);
save('water_ux_phi_.001_10hz.mat','ux');
save('water_uz_phi_.001_10hz.mat','uz');
save('water_p_phi_.001_10hz.mat','p');
save('water_uxt_phi_.001_10hz.mat','uxt');
save('water_uzt_phi_.001_10hz.mat','uzt');
save('water_pt_phi_.001_10hz.mat','pt');
save('water_recx_phi_.001_10hz.mat','recx');
save('water_recz_phi_.001_10hz.mat','recz');
%% gas
Kf(:)=9.6*10^7;
rhof(:)=600;
eta(:)=10^-5;
[ux,uz,p,uxt,uzt,pt,recx,recz,recp]= ... 
    poroelastic_solver_PLM(dims,C,Ks,rhos,Kf,rhof,eta,phi,k1,k3,R1,R3,rhoa,ome,lp,d0,vb);
save('gas_ux_phi_.001_10hz.mat','ux');
save('gas_uz_phi_.001_10hz.mat','uz');
save('gas_p_phi_.001_10hz.mat','p');
save('gas_uxt_phi_.001_10hz.mat','uxt');
save('gas_uzt_phi_.001_10hz.mat','uzt');
save('gas_pt_phi_.001_10hz.mat','pt');
save('gas_recx_phi_.001_10hz.mat','recx');
save('gas_recz_phi_.001_10hz.mat','recz');
%}
%{
figure(6)
set(gcf,'position',[80,80,1000,700]);
subplot(2,2,1)
imshow(real(oil_recz),.01*[min(real(oil_recz(:))),max(real(oil_recz(:)))]);
colorbar;
axis on;
%grid on;
ylim([0,dims.nt]);
xlabel({['x*' num2str(dims.dh) '[m]']});
ylabel({['t*' num2str(dims.dt) '[s]']});
ylim([0,250]);
title({['uzt']});

subplot(2,2,2)
imshow(real(water_recz),.01*[min(real(water_recz(:))),max(real(water_recz(:)))]);
colorbar;
axis on;
%grid on;
ylim([0,dims.nt]);
xlabel({['x*' num2str(dims.dh) '[m]']});
ylabel({['t*' num2str(dims.dt) '[s]']});
ylim([0,250]);
title({['uzt']});

subplot(2,2,3)
imshow(real(water_recz-oil_recz),.1*[min(real(water_recz(:)-oil_recz(:))),max(real(water_recz(:)-oil_recz(:)))]);
colorbar;
axis on;
%grid on;
ylim([0,dims.nt]);
xlabel({['x*' num2str(dims.dh) '[m]']});
ylabel({['t*' num2str(dims.dt) '[s]']});
ylim([0,250]);
title({['difference uzt']});

%print(gcf,['C:\Users\Yi\OneDrive\master thesis\VTI_poroelasticity\oil_10hz\rec_compare.png'],'-dpng','-r1200');
%}
