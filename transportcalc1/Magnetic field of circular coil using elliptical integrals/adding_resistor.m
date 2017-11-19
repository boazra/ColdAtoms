Itot=200;
Ir=0.05*Itot;
Iratio=(Itot-Ir)/Ir;
Rcoil=8.43e-3; %maybe smaller
R=Rcoil*Iratio;
Watt=R*Ir^2;

