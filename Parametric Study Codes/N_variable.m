% Before running the program, please import the WaterSat.xlsx file.
% Import as table, without renaming.
% Please check the end portion of the code to get desired result.

tci=20;
thi=50;
mh=0.6;
mc=0.6;
s=0.0035;
A=0.02306;
Ac=349.96*10^(-6);
Dh=6.756*10^(-3);
Dp=0.0254;
l=0.2159;
N=linspace(5,45,41);
Rf=0.00035;
t=0.0008;
km=237;

Q=size(40);
Tho=size(40);
Tco=size(40);
DelPh=size(40);
DelPc=size(40);

thoas=38;
tcoas=23;

for j=1:1:41
while 1 
    tavgh=(thi+thoas)/2;

for i=1:1:24
    if tavgh==WaterSat{i,1}
        dh=WaterSat{i,3};
        cph=WaterSat{i,4};
        kh=WaterSat{i,5};
        miuh=WaterSat{i,6}*(10^-3);
        prh=WaterSat{i,7};
        break
    elseif tavgh>WaterSat{i,1} && tavgh<WaterSat{i+1,1}  
        dh1=WaterSat{i,3};
        cph1=WaterSat{i,4};
        kh1=WaterSat{i,5};
        miuh1=WaterSat{i,6}*(10^-3);
        prh1=WaterSat{i,7};

        dh2=WaterSat{i+1,3};
        cph2=WaterSat{i+1,4};
        kh2=WaterSat{i+1,5};
        miuh2=WaterSat{i+1,6}*(10^-3);
        prh2=WaterSat{i+1,7};

        th1=WaterSat{i,1};
        th2=WaterSat{i+1,1};

        dh=((dh1-dh2)*(tavgh-th2))/(th1-th2)+dh2;
        cph=((cph1-cph2)*(tavgh-th2))/(th1-th2)+cph2;
        kh=((kh1-kh2)*(tavgh-th2))/(th1-th2)+kh2;
        miuh=((miuh1-miuh2)*(tavgh-th2))/(th1-th2)+miuh2;
        prh=((prh1-prh2)*(tavgh-th2))/(th1-th2)+prh2;
        break
    end
end

if rem(N(j),2)==0
        vh=(2*mh)/(dh*Ac*(N(j)+2));
else
        vh=(2*mh)/(dh*Ac*(N(j)+1));
end

Reh=(dh*Dh*vh)/miuh;

if Reh<100
        Gzh=(Reh*Prh)/(l/Dh);
        Nuh=1.86*(Gzh)^(1/3);
else
        Nuh=0.374*(Reh^0.668)*(prh^(1/3));
end

hh=(Nuh*kh)/Dh;

tavgc=(tci+tcoas)/2;

for i=1:1:24
    if tavgc==WaterSat{i,1}
        dc=WaterSat{i,3};
        cpc=WaterSat{i,4};
        kc=WaterSat{i,5};
        miuc=WaterSat{i,6}*(10^-3);
        prc=WaterSat{i,7};
        break
    elseif tavgc>WaterSat{i,1} && tavgc<WaterSat{i+1,1}  
        dc1=WaterSat{i,3};
        cpc1=WaterSat{i,4};
        kc1=WaterSat{i,5};
        miuc1=WaterSat{i,6}*(10^-3);
        prc1=WaterSat{i,7};

        dc2=WaterSat{i+1,3};
        cpc2=WaterSat{i+1,4};
        kc2=WaterSat{i+1,5};
        miuc2=WaterSat{i+1,6}*(10^-3);
        prc2=WaterSat{i+1,7};

        tc1=WaterSat{i,1};
        tc2=WaterSat{i+1,1};

        dc=((dc1-dc2)*(tavgc-tc2))/(tc1-tc2)+dc2;
        cpc=((cpc1-cpc2)*(tavgc-tc2))/(tc1-tc2)+cpc2;
        kc=((kc1-kc2)*(tavgc-tc2))/(tc1-tc2)+kc2;
        miuc=((miuc1-miuc2)*(tavgc-tc2))/(tc1-tc2)+miuc2;
        prc=((prc1-prc2)*(tavgc-tc2))/(tc1-tc2)+prc2;
        break
    end
end

if rem(N(j),2)==0
        vc=(2*mc)/(dc*Ac*(N(j)+2));
else
        vc=(2*mc)/(dc*Ac*(N(j)+1));
end

Rec=(dc*Dh*vc)/miuc;

if Rec<100
        Gzc=(Rec*Prh)/(l/Dh);
        Nuc=1.86*(Gzc)^(1/3);
else
        Nuc=0.374*(Rec^0.668)*(prc^(1/3));
end

hc=(Nuc*kc)/Dh;

U0=(1/hc+1/hh+2*Rf+t/km)^-1;

tmax=thi-tci;
cc=mc*cpc;
ch=mh*cph;
if cc>ch
    cmax=cc;
    cmin=ch;
else
    cmin=cc;
    cmax=cc;
end
qmax=cmin*tmax;
c=cmin/cmax;

entu=(U0*A*N(j))/cmin;

syms ntu0 entu0
ntu1=0.4;
er=0.001;
f=(ntu0-0.0166*(ntu0^2))-entu0;
g=diff(f,ntu0);
df_f=double(subs(g,{ntu0,entu0},{ntu1,entu}));
if df_f==0
    disp('NTU can not be determined');
else 
     while 1
        p=double(subs(f,{ntu0,entu0},{ntu1,entu}));
        if p==0
            break
        end
        df_f=double(subs(g,{ntu0,entu0},{ntu1,entu}));
        ntu2=ntu1-p/df_f;
        if abs((ntu2-ntu1)/ntu2)<=er
            ntu1=ntu2;
            break
        end
        ntu1=ntu2;
     end
end

ntu=ntu1;
if c==1
    e=ntu/(1+ntu);
else
    e=(1-exp(-ntu*(1-c)))/(1-c*exp(-ntu*(1-c)));
end
    
q=e*qmax;
tco=tci+q/cc;
tho=thi-q/ch;

if abs((tho-thoas)/thoas)<0.0001 && abs((tco-tcoas)/tcoas)<0.0001
    Tho(j)=tho;
    Tco(j)=tco;
    Q(j)=q;
    break
end

thoas=tho;
tcoas=tco;

end

if Reh>1 && Reh<10
    f=280/Reh;
elseif Reh>10 && Reh<100
    f=100/(Reh^0.589);
else
    f=12/(Reh^0.183);
end

vph=mh/(dh*((pi*Dp*Dp)/4));

DelPh(j)=(f*l/Dh)*(0.5*dh*vh*vh)+1.3*(0.5*dh*vph*vph);

if Rec>1 && Rec<10
    f=280/Rec;
elseif Rec>10 && Rec<100
    f=100/(Rec^0.589);
else
    f=12/(Rec^0.183);
end

vpc=mc/(dc*((pi*Dp*Dp)/4));

DelPc(j)=(f*l/Dh)*(0.5*dc*vc*vc)+1.3*(0.5*dc*vpc*vpc);

end

% disp('The values of Tho are');
% disp(Tho);
% disp('The values of Tco are');
% disp(Tco);
% disp('The values of Q are');
% disp(Q);
disp('The values of Hot Fluid Pressure Drops are');
disp(DelPh);
disp('The values of Cold Fluid Pressure Drops are');
disp(DelPc);

% plot(N,Tco,'-*',N,Tho,'-o');
% grid on;
% plot(N,Q,'-*');
% grid on
plot(N,DelPc,'-*',N,DelPh,'-o');
grid on


