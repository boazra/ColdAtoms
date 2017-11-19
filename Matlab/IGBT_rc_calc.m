%IGBT resistor and capacitor calculations

n = 3;
Vs_max = 8;
Is = 3e-3;
R_L = 1;
Vdc = 15;

Rs = (n*Vs_max-Vdc)/(n-1)/Is
Pd = Vs_max^2/Rs


