# OCT_klinearization_dispersion
Calibrating interpolation and dispersion compensation for OCT
How to properly perform k-linearization and dispersion compensation for OCT
1.Image two sets of Alines using a single reflector (such as the top surface of a glass slide) at two different depths, one close to the zero-delay position, the other around the middle of imaging depth.

2.Properly remove DC.
datB2 = zeros(size(datB)) ;
for ii = 1:Dim.nx
    datB2(:,ii) = double(datB(:,ii)) -smooth(double(datB(:,ii)),31);
end
The smooth function in MATLAB allows to remove the high frequency component, 31 is an empirial choice.

3.Check the phase of raw spectrum using hilbert transformation, select spectrum range with smooth phase accumulation.

In this example, the samples within [126,775] has smooth phase, so trim the samples outside of this range. Do the same for second set of Alines.

Take a look at the phase again.

There is an abrupt change around 100 pixels in the file2 phase, this is due to the low intensity in the spectrum where two SLDs overlap.
Without k-linearization and dispersion compensation, the FFT signal with have broadened PSFs and is depth dependent.

4.Subtract phaseB of file2 Alines by phaseA of file1 Alines. 

If there is abrupt change in either phase, select some Aline without abrupt change. Here the blue phase was having an abrupt change around 100 pixels, so I select an Aline from file2 without this change and do more smooth. Yellow line is the difference of these two phases.
5.Generate interpolation indices from the phase difference and interpolate both sets of Alines
%% subtract phaseA from phaseB to remove dispersion
dphase = smooth(phaseB-phaseA,11);
linPhase=linspace(dphase(1),dphase(end),length(dphase));
indice=find_interp_indice(dphase,linPhase);

%% interpolate file1 and file2 spectrum
datA5 = zeros(size(datA2)) ;
for ii = 1:Dim.nx
    datA5(:,ii)=get_interp(dphase,linPhase,datA2(:,ii),indice);
end
datB5 = zeros(size(datB2)) ;
for ii = 1:Dim.nx
    datB5(:,ii)=get_interp(dphase,linPhase,datB2(:,ii),indice);
end


After interpolation, PSFs become narrower, and depth dependency is removed.
6.Calculate the residual nonlinear phase, which is from dispersion.

The residual phase should be same for all depth. Here there is a difference, which is due to the abrupt change in second set of Alines.
7.Dispersion compensation through complex array.
datA7 = zeros(size(datA5)) ;
for ii = 1:Dim.nx
    datA7(:,ii)=datA5(:,ii).*exp(1j.*dphase1); % use dphase1 because dphase2 has abrupt changes
end
Here I used dphase1 because dphase2 has that abrupt change.
8.FFT of the final complex spectrum.

