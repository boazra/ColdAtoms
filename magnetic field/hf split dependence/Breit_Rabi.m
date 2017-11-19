gJ = 2.00233113;
gI = -0.0009951414;
h = 6.62606876 *10^(-34);
Ahfs = 3.41734130545215* 10^9;
Ehfs = Ahfs*2 ;%(*Hz*)
MuB = 1.399624624* 10^6;%(* Hz/G *)
Inuc = 3/2;
B1 = Ehfs/MuB/(gJ - gI);
x = (0:0.001:600 )/B1;
energy = zeros(2,5,size(x,2));

for F = Inuc-1/2:1:Inuc+1/2    
    for M = -F:F
        if F == Inuc - 1/2
           ene =  -(1/(2* (2* Inuc + 1))) - 1/2* sqrt(1 + (4* M)/(2 *Inuc + 1) *x + x.^2);
        else
            if M == -(Inuc + 1/2)
              ene = -(1/(2* (2* Inuc + 1))) + 1/2 *(1 - x);
            else
               ene = -(1/(2 *(2* Inuc + 1))) + 1/2* sqrt(1 + (4* M)/(2* Inuc + 1)* x + x.^2);
            end
        end 
        ene = ene +gI * M*x/(gJ - gI);
        energy(F,M+3,:) = ene * Ehfs/10^9;
    end
end

f1 = squeeze(energy(1,2:4,:))';
f2 = squeeze(energy(2,:,:))';

diffs = zeros(3,5,size(f1,1));
for i = 1:3
    for j = 1:5
        diffs(i,j,:) = f1(:,i)-f2(:,j);
    end
end


figure(1)
clf
hold on
plot(x*B1,f1,'r',x*B1,f2,'b')
xlabel('Magnetic Field (Gauss)')
ylabel('Energy')