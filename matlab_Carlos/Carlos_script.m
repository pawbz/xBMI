function Carlos_script
clc;
addpath(genpath('./'));
%%
a=colormap('gray');b=a(length(a):-1:1,:);groy=colormap(b);
P.perc=98;
%%
Diff=[250,500,750;450,450,500];ITRCLx=-100;ITRCLz=0;zSRC=250;
nx0=101;dx=10;
xRCV=(0:nx0-1).*dx;
xSRC=linspace(0,1000,200);    nsrc=length(xSRC);
dt=5e-4;
rcvz=10;rcvzz=repmat(rcvz,1,nx0);

nx=nx0;
rcvxtrue=xRCV;
rcvx=xRCV;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load('Cp.mat');   load('Ro.mat');
dz=5;  
I=M.*Ro;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mh=M(1,1).*ones(size(M));
Roh=Ro(1,1).*ones(size(M));
    [Lz,Lx]=size(M);itrcl=0;Lx=Lx-1;Lz=Lz-1;
    ilx=1+itrcl/dz;ilz=1+000/dz;
    m=Mh(ilz:ilz+Lz,ilx:ilx+Lx);ro1=Roh(ilz:ilz+Lz,ilx:ilx+Lx);
    x=ITRCLx-dz+(ilx:Lx+ilx).*dz;
    z=ITRCLz-dz+(ilz:Lz+ilz).*dz;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    modrezamp=2;iz1=1;
    vel=m(iz1:modrezamp:end,1:modrezamp:end);
    Rovel=ro1(iz1:modrezamp:end,1:modrezamp:end);
    x=x(1:modrezamp:end);z=z(iz1:modrezamp:end);dz2=dz*modrezamp;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Y=diff(I);Yy=diff(I.').';
    I=(abs(Y(:,1:end-1)+Yy(1:end-1,:)));

    h1=figure(90);
    set(h1,'position',[50 350 400 200]);
       figure(90);clf;imagesc(x,z,I);colorbar;
           hold on;plot(xRCV(1:1:end),rcvz,'kv','Markerfacecolor','y','markersize',12);
           hold on;plot(xSRC(1:1:end),zSRC(1:1:end),'kp','Markerfacecolor','r','markersize',5);
           caxis([1500,6000]);
           cok=colormap('gray');       colormap(cok(end:-1:1,:))
           axis equal;axis tight;
           colorbar;
           xlabel('Position (m)','fontname','timesnewroman');
           ylabel('Depth (m)'   ,'fontname','timesnewroman');
           hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist('vzc','var')==0,load('vz_t.mat');end
nt=size(vzc,1);

ttaxis=(-nt:nt-1).*dt;

dom=2*pi/(nt*dt);nf=ceil(nt/2);om=(0:nf-1)*dom;f=om/2/pi;df=dom/2/pi;

fdnf=1000;
dnf=find(abs(f-fdnf)==min(abs(f-fdnf)));

Nsrc=1:nsrc;

    f2=linspace(0,1/dt/2,nt);
    Taper=repmat(Gibbbs(f2,5,15,900,1000),1,nx0);

    Ccc=zeros(nx0,nx0,nt);
    Vzc=zeros(nt,nx0,length(Nsrc));
    c=zeros(2*nt,nx0,nx0);
    is1=0;
    for is=Nsrc
        is1=is1+1;
        Vzc(:,:,is1)=ftfour([vzc(:,:,is1);zeros(nt,nx0)],nt);
    end
    Vzc=permute(Vzc,[2,3,1]);
    for iw=1:dnf
        Vu=[Vzc(:,:,iw),Vzc(:,:,iw)];
        Ccc(:,:,iw)=(Vu*Vu').';
    end
    Ccc2=permute(Ccc,[3,1,2]);
    for ivs=1:nx0
        c(:,:,ivs)=fftshift(itfour(Taper.*Ccc2(:,:,ivs),2*nt),1);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ivs=round(nx/2);
%     for ivs=1:nx0
        tmp=squeeze(c(:,ivs,:));
        figure(232);
        imagesc(rcvxtrue,ttaxis,tmp);
        clip = prctile(abs(tmp(:)),[P.perc]);caxis([-clip, clip]);
        colormap(groy);
        hold on;
        plot(rcvxtrue(ivs),0,'pk','markerfacecolor','r');
        hold off;
        axis([rcvxtrue(1) rcvxtrue(end) -0.6 0.9]);pause(0.0001);
%     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotmigr=1;
dtt=dt;
%% Wavelet %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
FwF=zeros(2*nt,1);FwF(nt+1)=1;
%% Model
dz=dz2;
v=vel;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nrcvB=find(abs(rcvx-0)==min(abs(rcvx-0)));
nrcvE=find(abs(rcvx-1000)==min(abs(rcvx-1000)));
vsrcB=find(abs(rcvx-0)==min(abs(rcvx-0)));
vsrcE=find(abs(rcvx-1000)==min(abs(rcvx-1000)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt=0.0005;     irsample=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   rcv_s=nrcvB:1:nrcvE;%[2:2:82];
%    rcv_s=76;
   rcv=rcv_s;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   RcvTaper=Gibbbs(rcv_s,rcv_s(1)-10,rcv_s(1)+20,rcv_s(end)-20,rcv_s(end)+10);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     vsrc=vsrcB:4:vsrcE;
    vsrc=51;
    vsource=rcvx(vsrc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Frequency values for migration
f1=2;
f2f=115;
f2=65;
iidifb=2;
iidiff=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dtime=0.005;
%% Elliptic migration times
% time=-0.0:dtime:1.4;
%% Hyperbollic migration times
time=-0.45:dtime:0.1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_i=round(time./dt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    t_b_coor=zeros(length(rcv),1);
    cut=length(v)-size(v,1);    vzeros=[v;zeros(cut,length(v))]./(max(x)-min(x));
    T1bm=zeros(size(v(:),1),length(rcv));eTTbm=T1bm;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i1r=1:length(rcv)
        ir=rcv(i1r); fprintf(1,'\nrcv = %d\n',ir);            
        ircvx=find(abs(x-rcvx(ir))==min(abs(x-rcvx(ir))),1);            
        ircvz=find(abs(z-rcvzz(ir))==min(abs(z-rcvzz(ir))),1);
        v0=vel(ircvz,ircvx);

        boundary = [ircvz,ircvx,0];
        t_temp   = eikonal(vzeros,boundary);
        trcv=(t_temp(1:size(vel,1),:)).*1;
        T1b=(-1i./(4.*pi.*v0).*sqrt(v)).*(v.*v0.*sqrt(-0.5.*1j.*Rovel));
        eTTb=exp(trcv);
        t_b_coor(i1r)=round(max(trcv(:))/dt);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        T1bm(:,i1r)=single(T1b(:));
        eTTbm(:,i1r)=single(eTTb(:));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      

for i1vsrc=1:length(vsrc)
    ivsrc=vsrc(i1vsrc);
        fprintf(1,'virtual source at = %d m\n',vsource(i1vsrc)); 
        isrc=find(rcvx==vsource(i1vsrc));
        srcx=vsource(i1vsrc);
        srcz=rcvzz(isrc);
        v0=vel(find(abs(z-srcz)==min(abs(z-srcz)),1,'first'),find(abs(x-srcx)==min(abs(x-srcx)),1,'first'));
        
        ivsrcx=find(abs(x-rcvx(isrc))==min(abs(x-rcvx(isrc))),1);            
        ivsrcz=find(abs(z-rcvzz(isrc))==min(abs(z-rcvzz(isrc))),1);
        fprintf(1,'source beam iso \n'); 
        boundary = [ivsrcz,ivsrcx,0];
        t_temp   = eikonal(vzeros,boundary);
        trcv=(t_temp(1:size(vel,1),:));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        T1f=(-1i./(4.*pi.*v0).*sqrt(v)).*(v.*v0.*sqrt(-0.5.*1j.*Rovel));
        eTTf=exp(trcv);
        t_f_coor=round(max(trcv(:))/dt);
        t_f=max(trcv(:));
        T1fm=single(T1f(:));
        eTTfm=single(eTTf(:));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        uI=zeros(size(v(:)));uId=uI;
        ttaxis_resamp=(-nt:1*nt-1).*dtt;
        ttaxis_resamp1=resample( ttaxis_resamp ,irsample,1);
        taxis_resamp1=ttaxis_resamp1(nt+1:2*nt);
        f_resamp=       (resample( FwF ,irsample,1));
        c_resamp=squeeze(resample( c(:,ivsrc,:) ,irsample,1));
        itmax=max([t_f_coor;t_b_coor(:)]);
        RcvTaper2=repmat(RcvTaper,1,itmax);size(RcvTaper2)

        C_f=zeros(1*itmax,1);
        C_b=zeros(length(rcv),itmax);
        U_cb=Gibbbs(1:itmax,21,51,itmax,itmax+1);
        newfmax=1/2/(1*dtt/irsample);
        om = 2*pi*linspace(0,newfmax,1*itmax);
        om1=om;om1(1)=om1(2)./2;
        Invomb=(1./(-1j.*om1));     
        Invomf=ones(nt,1);
        U=Gibbbs(om/2/pi,0.5,8,newfmax,newfmax+1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for it=1:length(t_i)
        pause(0.001);
        tic
        fprintf(2,'\ntime = %.2f\n',time(it));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        shifter_for=t_i(it);
        shifter_back=t_i(it);   

        t_coor  =t_f_coor;
        fill_zeros=itmax-t_coor;
        %% Hyper
        temp=circshift(f_resamp,-shifter_for);
        itf=nt+(1:t_coor);
        %% Elipt
    %         temp=circshift((f_resamp),shifter_for); %FwF starts at xs at 9s
    %         itf=(1.*nt+(1:t_coor));
        C_f=[(temp(itf));zeros(fill_zeros,1)];

        for i1r=1:length(rcv)
            ir=rcv((i1r));
            t_coor=max(t_b_coor(i1r,:));
            fill_zeros=itmax-t_coor;
            temp=circshift(c_resamp,-shifter_back);
            itb=nt+(1:t_coor);
            C_b(i1r,:)=U_cb.*[(temp(itb,ir));zeros(fill_zeros,1)];
        end

        greenf1 = fft(C_f,(2*itmax),1);greenf=(greenf1(1:1*itmax,:).*U);
        greenb1 = fft(C_b,(2*itmax),2);greenb=RcvTaper2.*greenb1(:,1:1*itmax).*repmat(U.',[length(rcv) 1]);

        greenb=single(greenb);
        om=single(om);
        Invomb=single(Invomb);
        Invomf=single(Invomf);
        if1=find(abs(2*pi*f1-om)==min(abs(2*pi*f1-om)));
        if2=find(abs(2*pi*f2-om)==min(abs(2*pi*f2-om)));
        if2f=find(abs(2*pi*f2f-om)==min(abs(2*pi*f2f-om)));
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FREQUENCY Operatia %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        axF=[0 1000 0 700];
        h1=figure(23332);
        if isempty(findobj('Name','Pawan'))
            set(h1,'Position',[10 50 1200 300],'Name','Pawan');
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%
        Nf=if1:iidifb:if2;
        uub=zeros(size(v(:)));
        urb=zeros(size(v(:),1),length(rcv));
        for iw=Nf
            clear urbw
            func=repmat(squeeze((greenb(:,iw))),1,size(T1bm,1)).';
            urbw=Invomb(iw).*(T1bm).*(eTTbm).^(1j*om(iw)).*func;urbw(isnan(urbw))=0;
            urb=urbw+urb;
        end
        u1b=sum(urb(:,:),2);
        uub=uub+u1b;
        if plotmigr==1
            figure(h1);subplot(1,3,2);
            imagesc(x,z,real(reshape(uub,size(vel))));
            tmp = prctile(abs(uub(:)),[P.perc]);caxis([-tmp, tmp]);
            colormap(groy);axis equal;axis tight;
            hold on;plot(Diff(1,:),Diff(2,:),'b+');hold off
            axis(axF);
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%
        Nf=if1:iidiff:if2f;
        urf=zeros(size(v(:)));
        urfa=zeros(size(v(:),1),1);
        for iw=Nf
            clear urfw
            urfw=Invomf(iw).*T1fm.*eTTfm.^(1j.*om(iw)).*(greenf(iw));urfw(isnan(urfw))=0;
            urfa=(urfw)+urfa;
        end
        urf=urf+urfa;
        if plotmigr==1
            figure(h1);subplot(1,3,1);
            imagesc(x,z,real(reshape(urf,size(vel))));
            tmp = prctile(abs(urf(:)),[P.perc]);caxis([-tmp, tmp]);
            colormap(groy);axis equal;axis tight;axis equal;axis tight;
            hold on;plot(Diff(1,:),Diff(2,:),'b+');hold off
            axis(axF);
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%
        uI2=urf.*uub;
        uI=uI+uI2;
%%%%%%%%%%%%%%%%%%%%%%%%%%
        if plotmigr==1
            figure(h1);subplot(1,3,3);
            imagesc(x,z,real(reshape(uI,size(vel))))
            tmp = prctile(abs(uI(:)),[P.perc]);caxis([-tmp, tmp]);
            colormap(groy);axis equal;axis tight;
            hold on;plot(Diff(1,:),Diff(2,:),'b+');hold off
            axis(axF);
        end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        toc
    end
end    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

function [u]=Gibbbs(t,tol,tcl,tc,to)
    u=zeros(length(t),1);
    for ik=1:length(t)
        if t(ik)<tol || t(ik)>to
            u(ik)=0;
        elseif t(ik)<tcl && t(ik)>=tol
            u(ik)=0.5+0.5*cos(pi*(1/abs(tol-tcl))*(t(ik)-tcl));
        elseif t(ik)<tc && t(ik)>=tcl
            u(ik)=1;
        elseif t(ik)>=tc && t(ik)<=to
            u(ik)=0.5+0.5*cos(pi*(1/abs(tc-to))*(t(ik)-tc));
        end
    end
end

function [U]=ftfour(u,nfc)

    [nf,nrcv]=size(u);
    U=zeros(nfc,nrcv);

    for ir=1:nrcv
        tmp=fft(u(:,ir),nf);
        U(:,ir)=tmp(1:nfc);
    end
    clear temp;
end

function [u]=itfour(U,nsp)

    [nfc,nrcv] = size(U);
    tmp=zeros(nsp,1);
    u=zeros(nsp,nrcv);

    for ir=1:nrcv
        tmp(1:nfc)=U(:,ir);             
        tmp(nsp-nfc+2:nsp) = conj(U(nfc:-1:2,ir));
        u(:,ir) = real(ifft(tmp));
    end
    clear tmp;
end