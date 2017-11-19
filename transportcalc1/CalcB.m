function [Bx,By,Bz,B_mag] = CalcB(coil,sz,sx)
%calculates the magnetic field of 2 coils in anti helmholtz configuration
%with 1 amper current. in the area of interest sz,sx, coils are at coil.x0
%horizontally.
%   Detailed explanation goes here
%%%%% Quad coils
%%
%addpath('/ph2users/hagaie/Documents/MATLAB/Modified Transport/new calc/Magnetic field of circular coil using elliptical integrals')
global u0
u0=4*pi*1e-7;
 %new quad coils
% coil.Sr=1e-3; %radial spacing of windings
% coil.Sz=2.6e-3; %axial spacing of windings
% coil.Rin=15e-3; %inner radius
% coil.Zmin=18e-3; %lower coil axial distance 14
% coil.Nz=4; %windings in axial direction - 4
% coil.Nr=25; %windings in radial direction - 6
Iq=1; %Quad current - Ampere

% % Transport coils
% Sr=1e-3; %radial spacing of windings
% Sz=2.5e-3; %axial spacing of windings
% Rin=9e-3; %inner radius
% Zmin=16.5e-3; %lower coil axial distance 14
% Nz=2; %windings in axial direction - 4
% Nr=17; %windings in radial direction - 6
% Iq=-50; %Quad current - Ampere


% sz=(-130:1:107)*1e-3;
% sx=(-129:1:129)*1e-3;
[X,Z]=meshgrid(sx,sz);
%axes- Z - quad symmetry +gravity, Y - Ioffe symmtery
Bx=zeros(size(X));By=Bx;Bz=Bx;

for j=1:coil.Nr
    for k=1:coil.Nz
        [Bx1,By1,Bz1]=magnetic_field_current_loop(X,0,Z,coil.x0,0,coil.Zmin+(k-0.5)*coil.Sz,coil.Rin+(j-0.5)*coil.Sr,Iq);
        Bx=Bx+Bx1;
        By=By+By1;
        Bz=Bz+Bz1;
        [Bx1,By1,Bz1]=magnetic_field_current_loop(X,0,Z,coil.x0,0,-coil.Zmin-(k-0.5)*coil.Sz,coil.Rin+(j-0.5)*coil.Sr,-Iq);
        Bx=Bx+Bx1;
        By=By+By1;
        Bz=Bz+Bz1;
    end
end
B_mag=sqrt(Bx.^2+By.^2+Bz.^2);
% %% plotting
% figure(24)
% imagesc(sx,sz,B_mag*1e4)
% axis image
% % figure(5)
% % plot(sx,B_magQ(281,:)*1e4)



end

