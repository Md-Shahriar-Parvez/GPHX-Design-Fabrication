% Before running the program, please import the WaterSat.xlsx file.
% Import as table, without renaming.

% Input
tci=20;
thi=50;
tho=45.55;
mh=0.6;
mc=0.6;
% Geometric Parameters and Properties
s=0.0035;
A=0.02306;
Ac=349.96*10^(-6);
Dh=6.756*10^(-3);
Dp=0.0254;
l=0.2159;
Rf=0.00009;
t=0.0008;
km=237;
% Thermophysical Properties of Hot Fluid
tavgh=(thi+tho)/2;
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
% Thermophysical Properties of Cold Fluid
q=mh*cph*(thi-tho);
tcoas=25;
tavgc=(tci+tcoas)/2;
for i=1:1:24
    if tavgc==WaterSat{i,1}
        cpc=WaterSat{i,4};
        break
    elseif tavgc>WaterSat{i,1} && tavgc<WaterSat{i+1,1}  
        cpc1=WaterSat{i,4};
        cpc2=WaterSat{i+1,4};       
        tc1=WaterSat{i,1};
        tc2=WaterSat{i+1,1};
        cpc=((cpc1-cpc2)*(tavgc-tc2))/(tc1-tc2)+cpc2;
        break
    end
end       
    q0=mc*cpc*(tcoas-tci);
    if q0==q
        tco=tcoas;
    else
        itr1=0;
        while 1
            tcoas=q/(mc*cpc)+tci;
            tavgc=(tci+tcoas)/2;
            for i=1:1:24
    if tavgc==WaterSat{i,1}
        cpc=WaterSat{i,4};
        break
    elseif tavgc>WaterSat{i,1} && tavgc<WaterSat{i+1,1}  
        cpc1=WaterSat{i,4};
        cpc2=WaterSat{i+1,4};      
        tc1=WaterSat{i,1};
        tc2=WaterSat{i+1,1};
        cpc=((cpc1-cpc2)*(tavgc-tc2))/(tc1-tc2)+cpc2;
        break
    end
            end
            q0=mc*cpc*(tcoas-tci);
            itr1=itr1+1;
            if abs(q0-q)>0 && abs(q0-q)<0.5
            tco=tcoas;
            break
            end
        end
    end 
    tavgc=(tci+tco)/2;
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
% Calculation of Required Number of Plates and Thermal Performance
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
e=(q/qmax)*100;
c=cmin/cmax;
syms ntu0 c0 e0
ntu1=0.25;
e=e/100;
er=0.001;
itr2=0;
if c==1
    f=ntu0/(1+ntu0)-e0;
else
    f=(1-exp(-ntu0*(1-c0)))/(1-c0*exp(-ntu0*(1-c0)))-e0;
end
g=diff(f,ntu0);
df_f=double(subs(g,{ntu0,c0,e0},{ntu1,c,e}));
if df_f==0
    disp('NTU can not be determined');
else 
    while 1
        itr2=itr2+1;
        p=double(subs(f,{ntu0,c0,e0},{ntu1,c,e}));
        if p==0
            break
        end
        df_f=double(subs(g,{ntu0,c0,e0},{ntu1,c,e}));
        ntu2=ntu1-p/df_f;
        if abs((ntu2-ntu1)/ntu2)<=er
            ntu1=ntu2;
            break
        end
        ntu1=ntu2;
    end
end
fcorr=1-0.0166*ntu1;
ntu1=fcorr*ntu1;
UN=(ntu1*cmin)/A;
N=12;
itr3=0;
while 1
    itr3=itr3+1;
    U=UN/N;
    N=ceil(N);
    if rem(N,2)==0
        vh=(2*mh)/(dh*Ac*(N+2));
        vc=(2*mc)/(dc*Ac*(N+2));
    else
        vh=(2*mh)/(dh*Ac*(N+1));
        vc=(2*mc)/(dc*Ac*(N+1));
    end
    Reh=(dh*Dh*vh)/miuh;
    Rec=(dc*Dh*vc)/miuc;
    if Reh<100
        Gzh=(Reh*Prh)/(l/Dh);
        Nuh=1.86*(Gzh)^(1/3);
    else
        Nuh=0.374*(Reh^0.668)*(prh^(1/3));
    end
    if Rec<100
        Gzc=(Rec*Prc)/(l/Dh);
        Nuc=1.86*(Gzc)^(1/3);
    else
        Nuc=0.374*(Rec^0.668)*(prc^(1/3));
    end
    hh=(Nuh*kh)/Dh;
    hc=(Nuc*kc)/Dh;
    U0=(1/hc+1/hh+2*Rf)^-1; 
    if abs((U0-U)/U)<=0.01
        Ns=N;
        U1=UN/Ns;
        break
    end
    N=UN/U0;
end
% Output
disp('All values are in SI unit');
fprintf('ρh=%.4f\n',dh);
fprintf('cph=%.4f\n',cph);
fprintf('kh=%.4f\n',kh);
fprintf('μh=%.4f\n',miuh);
fprintf('Prh=%.4f\n',prh);
fprintf('ρc=%.4f\n',dc);
fprintf('cpc=%.4f\n',cpc);
fprintf('kc=%.4f\n',kc);
fprintf('μc=%.4f\n',miuc);
fprintf('Prc=%.4f\n',prc);
fprintf('hh=%.4f\n',hh);
fprintf('hc=%.4f\n',hc);
fprintf('Uo=%.4f\n',U1);
fprintf('Q=%.4f\n',q);
fprintf('Tco=%.4f\n',tco);
fprintf('Required number of plates=%.4f\n',Ns);




