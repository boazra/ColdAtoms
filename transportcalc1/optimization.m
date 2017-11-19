%optimization

%%perferct quadrupole
sx=(-10:0.1:10)*1e-3;
sz=sx;
[X,Z]=meshgrid(sx,sz);
dBdz=0.1; %1 T/m =100 G/cm
dBdx=-dBdz/3;
Bx=dBdx*X;
Bz=dBdz*Z;
Btarget=sqrt(Bx.^2+Bz.^2);
clear 'Bz' 'Bx''X' 'Z'
%% real calc
%sz=(-10:0.1:10)*1e-3;
% sx=(340:0.025:360)*1e-3;

Quad.Sr=1e-3; %radial spacing of windings
Quad.Sz=2.5e-3; %axial spacing of windings
Quad.Rin=15e-3; %inner radius
Quad.Zmin=18e-3; %lower coil axial distance 14
Quad.Nz=4; %windings in axial direction - 4
Quad.Nr=25; %windings in radial direction - 6
Quad.x0=0.423;


TC2.Sr=1e-3; %radial spacing of windings
TC2.Sz=2.5e-3; %axial spacing of windings
TC2.Rin=9e-3; %inner radius
TC2.Zmin=28e-3; %lower coil axial distance 14
TC2.Nz=2; %windings in axial direction - 4
TC2.Nr=17; %windings in radial direction - 6
TC2.x0=0.38;

TC1.Sr=1e-3; %radial spacing of windings
TC1.Sz=2.5e-3; %axial spacing of windings
TC1.Rin=9e-3; %inner radius
TC1.Zmin=22e-3; %lower coil axial distance 14
TC1.Nz=2; %windings in axial direction - 4
TC1.Nr=17; %windings in radial direction - 6
TC1.x0=0.35;
Iqmax=4;
IT1max=12;
IT2max=12;
Ires=0.2;
Iq=0:Ires:Iqmax;nq=length(Iq);
IT1=0:Ires:IT1max;nT1=length(IT1);
IT2=0:Ires:IT2max;nT2=length(IT2);
 
x0=350e-3:10e-3:423e-3;
I1=zeros(length(x0),1);I2=I1;I3=I1;m=I1;
for l=1:length(x0)

sx=x0(l)+sz;
[BxQ,ByQ,BzQ]=CalcB(Quad,sz,sx);
[BxT2,ByT2,BzT2]=CalcB(TC2,sz,sx);
[BxT1,ByT1,BzT1]=CalcB(TC1,sz,sx);

OP=zeros(nT1,nT2,nq);
for i=1:nT1
    for j=1:nT2
        for k=1:nq
            B=sqrt((IT1(i)*BxT1+IT2(j)*BxT2+Iq(k)*BxQ).^2+(IT1(i)*ByT1+IT2(j)*ByT2 ...
                +Iq(k)*ByQ).^2+(IT1(i)*BzT1+IT2(j)*BzT2+Iq(k)*BzQ).^2);
            OP(i,j,k)=sum(sum((B-Btarget).^2));
        end
    end
end

[m(l),in]=min(OP(:));
[im,jm,km]=ind2sub(size(OP),in);
I1(l)=IT1(im);
I2(l)=IT2(jm);
I3(l)=Iq(km);
end
save('currents_12triplet_newcalc_dBdz10Gcm_AR2.mat','I1','I2','I3','m')
figure(1);clf;plot(x0,I1,x0,I2,x0,I3)
clock
%  Bf=sqrt((I1(l)*BxT1+I2(l)*BxT2+I3(l)*BxQ).^2+(I1(l)*ByT1+I2(l)*ByT2 ...
%                 +I3(l)*ByQ).^2+(I1(l)*BzT1+I2(l)*BzT2+I3(l)*BzQ).^2);
% figure(2)
% clf
% subplot(1,2,1)
% imagesc(sx,sz,Bf*1e4)
% subplot(1,2,2)
% imagesc(sx,sz,Btarget*1e4)
% figure(3)
% clf
% plot(sz,Bf(:,401)*1e4)
% hold on
% plot(sz,Bf(401,:)*1e4)
% hold off