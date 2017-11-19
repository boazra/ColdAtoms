
function AtomNumber = CountAtomsAbsorption(IMG,Idark)
    Isat = 1.75;
    I = 0.05*Isat;
    ODsat = 6.0;

    IMGd = double(IMG);
    Intensity_ratio = (IMGd(:,:,2)-Idark)./(IMGd(:,:,1)-Idark);
    Intensity_ratio(Intensity_ratio < 1 ) = 1;
    ODmeas = reallog(Intensity_ratio );
    ODmeas(isinf(ODmeas)) = 0;
    ODmeas(isnan(ODmeas)) = 0;
    % Remove noise by looking at difference between nearest neighbors
    ODmeas2 = Smooth(ODmeas,6);
    ODmeas = Smooth(ODmeas2,0.8);

    %calculate from Lewandowski
    temp = (1-exp(-ODsat))./(exp(-ODmeas)-exp(-ODsat));
    temp(temp<1) = 1;
    ODmod = log(temp);
    ODactual = ODmod + (1+exp(-ODmod))*I/Isat;
    lambda = 766e-9;
    detuning = 0;
    linewidth  =2*pi*6.03e6;
    branching_ratio = 1;
    A = (branching_ratio/2)*(3*lambda^2/2/pi)*(1+4*(detuning/linewidth)^2)^-1;
    AtomCount = ODactual/A;

    % extract only area with atoms
    total_atom_slice = sum(AtomCount);
    thresh = 0.8;
    peak_window = total_atom_slice(round(size(total_atom_slice,2)/4):round(size(total_atom_slice,2)*3/4));
    [maxval, ~] = max(peak_window);
    delta = peak2peak(peak_window);
    I = find(total_atom_slice<maxval-thresh*delta);
    di = diff(I);
    [~, arg] = max(di);
    interesting_window = AtomCount(:,I(arg):I(arg+1));
    pixel_size = 0.0055e-3^2;
    AtomCount = interesting_window * pixel_size;
    AtomNumber = sum(sum(AtomCount));
end
