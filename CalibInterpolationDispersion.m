%% Generate interpolation indices and dispersion phases for OCT using two sets of Alines of single reflector at different depths
%% Select file location %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
datapath  = 'D:\SDOCT_Klinear\';
cd(datapath)
% get data information                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
Dim.nx = 400;
Dim.ny = 1;
Dim.nk = 1024;
Dim.nxRpt = 1;
Dim.nyRpt = 1;
zRange = round(Dim.nk/2);
% Data processing info
files = dir('*.bin');
nFile = length(files);  % number of file

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% process first file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename=files(1).name;
filePath=[datapath,filename];
% read spectrum
fid=fopen(filePath,'r','l');
datA = fread(fid, Dim.nk*Dim.nxRpt*Dim.nx*Dim.nyRpt, 'uint16');
datA=reshape(datA, [Dim.nk Dim.nxRpt*Dim.nx Dim.nyRpt]);

figure;
for ii = 2:50:Dim.nx
    plot(datA(:,ii),'LineWidth',1);hold on
end
xlim([10,Dim.nk])
title('file1 raw spectrum');ax=gca;ax.FontSize=18;
%%%%%%%%%%%%%%%%%%%%%%%%%%% need to smooth the mean otherwise some
%%%%%%%%%%%%%%%%%%%%%%%%%%% signal will get subtracted out
datA2 = zeros(size(datA)) ;
for ii = 1:Dim.nx
    datA2(:,ii) = double(datA(:,ii)) -smooth(datA(:,ii),11);
end

figure;
for ii = 2:50:Dim.nx
    plot(datA2(:,ii),'LineWidth',1);hold on
end
xlim([1,Dim.nk])
title('file1 DC removed raw spectrum');ax=gca;ax.FontSize=18;

%% Display phase using hilbert transformation
datA4 = zeros(size(datA2)) ;
for ii = 1:Dim.nx
    tmp=hilbert(datA2(:,ii));
    datA4(:,ii)=unwrap(angle(tmp));
end
figure
for ii = 2:50:Dim.nx
    plot(datA4(:,ii),'LineWidth',1);hold on
end
ylabel('rad');xlim([1,Dim.nk])
title('file1 phase after HT');ax=gca;ax.FontSize=18;
%% trim non-interference samples
datA2=datA2(101:700,:);
Dim.nk = size(datA2,1);
zRange = round(Dim.nk/2);

figure
for ii = 2:50:Dim.nx
    plot(datA2(:,ii),'LineWidth',1);hold on
end
xlim([1,Dim.nk]);title('file1 trimmed spectrum');ax=gca;ax.FontSize=18;

RR0 = ifft(datA2,[],1);
RRA = abs(RR0(1:zRange,:,:));
figure
for ii = 2:50:Dim.nx
    plot(RRA(:,ii),'LineWidth',2);hold on
end
xlim([1,zRange]);title('file1 raw Alines');ax=gca;ax.FontSize=18;
%% re-look at phase after trimming spectrum
datA4 = zeros(size(datA2)) ;
for ii = 1:Dim.nx
    tmp=hilbert(datA2(:,ii));
    datA4(:,ii)=unwrap(angle(tmp));
end
figure
for ii = 2:50:Dim.nx
    plot(datA4(:,ii),'LineWidth',1);hold on
end
ylabel('rad');xlim([1,Dim.nk])
title('file1 phase after spectrum trim');ax=gca;ax.FontSize=18;
%% smooth phase1
phaseA = smooth(mean(datA4,2),21);figure;plot(phaseA,'LineWidth',2);title('phaseA');ax=gca;ax.FontSize=18;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% process second file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename=files(9).name;
filePath=[datapath,filename];
Dim.nk = 1024;
% read spectrum
fid=fopen(filePath,'r','l');
datB = fread(fid, Dim.nk*Dim.nxRpt*Dim.nx*Dim.nyRpt, 'uint16');
datB=reshape(datB, [Dim.nk Dim.nxRpt*Dim.nx Dim.nyRpt]);
figure;
for ii = 2:50:Dim.nx
    plot(datB(:,ii),'LineWidth',1);hold on
end
xlim([1,Dim.nk])
title('file2 raw spectrum');ax=gca;ax.FontSize=18;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% need to smooth the mean otherwise some
%%%%%%%%%%%%%%%%%%%%%%%%%%% signal will get subtracted out
datB2 = zeros(size(datB)) ;
for ii = 1:Dim.nx
    datB2(:,ii) = double(datB(:,ii)) -smooth(datB(:,ii),11);
end
% trim file2 spectrum to be the same as file1 spectrum
datB2=datB2(101:700,:);
Dim.nk = size(datB2,1);
figure
for ii = 2:50:Dim.nx
    plot(datB2(:,ii),'LineWidth',1);hold on
end
xlim([1,Dim.nk]);title('file2 trimmed spectrum');ax=gca;ax.FontSize=18;

RR0 = ifft(datB2,[],1);
RRB = abs(RR0(1:zRange,:,:));
figure
for ii = 2:50:Dim.nx
    plot(RRB(:,ii),'LineWidth',2);hold on
end
xlim([1,zRange]);title('file2 raw Alines');ax=gca;ax.FontSize=18;
%% Display phase using hilbert transformation
datB4 = zeros(size(datB2)) ;
for ii = 1:Dim.nx
    tmp=hilbert(datB2(:,ii));
    datB4(:,ii)=unwrap(angle(tmp));
end
figure
for ii = 2:50:Dim.nx
    plot(datB4(:,ii),'LineWidth',1);hold on
end
ylabel('rad');xlim([1,Dim.nk])
title('file2 phase after HT');ax=gca;ax.FontSize=18;

%% smooth phase2. carefully choose Aline without abrupt change in phase
phaseB = smooth(mean(datB4,2),21);figure;plot(phaseB,'LineWidth',2);title('phaseB');ax=gca;ax.FontSize=18;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% k linearization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;plot(phaseA,'red','LineWidth',2);hold on;plot(phaseB,'blue','LineWidth',2);
title('phase of two Alines at different depths');ax=gca;ax.FontSize=18;
% subtract phaseA from phaseB to remove dispersion
dphase = smooth(phaseB-phaseA,11);
linPhase=linspace(dphase(1),dphase(end),length(dphase));
figure;plot(dphase,'LineWidth',2);hold on;plot(linPhase,'LineWidth',2);title('phaseB');ax=gca;ax.FontSize=18;
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
figure
for ii = 2:50:Dim.nx
    plot(datA5(:,ii),'LineWidth',2);hold on
end
xlim([1,Dim.nk])
title('file1 spectrum after interpolation');ax=gca;ax.FontSize=18;
figure
for ii = 2:50:Dim.nx
    plot(datB5(:,ii),'LineWidth',2);hold on
end
xlim([1,Dim.nk])
title('file2 spectrum after interpolation');ax=gca;ax.FontSize=18;

RR0 = ifft(datA5,[],1);
RRA = abs(RR0(1:zRange,:,:));
figure
for ii = 2:50:Dim.nx
    plot(RRA(:,ii),'LineWidth',2);hold on
end
xlim([1,zRange]);title('file1 Alines after interpolation');ax=gca;ax.FontSize=18;
RR0 = ifft(datB5,[],1);
RRB = abs(RR0(1:zRange,:,:));
figure
for ii = 2:50:Dim.nx
    plot(RRB(:,ii),'LineWidth',2);hold on
end
xlim([1,zRange]);title('file2 Alines after interpolation');ax=gca;ax.FontSize=18;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% dispersion compensation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
datA6 = zeros(size(datA5)) ;
for ii = 1:Dim.nx
    tmp=hilbert(datA5(:,ii));
    datA6(:,ii)=unwrap(angle(tmp));
end
phase1 = smooth(mean(datA6,2));
line1 = linspace(phase1(1),phase1(end),length(phase1))';
dphase1 = phase1-line1;figure;plot(dphase1,'LineWidth',2);hold on;
datB6 = zeros(size(datB5)) ;
for ii = 1:Dim.nx
    tmp=hilbert(datB5(:,ii));
    datB6(:,ii)=unwrap(angle(tmp));
end
phase2 = smooth(mean(datB6,2));
line2 = linspace(phase2(1),phase2(end),length(phase2))';
dphase2 = phase2-line2;plot(dphase2,'LineWidth',2);hold on;
xlim([1,Dim.nk]);title('residual dispersion');ax=gca;ax.FontSize=18;
%% dispersion compensation and fft
datA7 = zeros(size(datA5)) ;
for ii = 1:Dim.nx
    datA7(:,ii)=datA5(:,ii).*exp(1j.*dphase1); % use dphase1 because dphase2 has abrupt changes
end
RR0 = ifft(datA7,[],1);
RRA = abs(RR0(1:zRange,:,:));
figure
for ii = 2:50:Dim.nx
    plot(RRA(:,ii),'LineWidth',2);hold on
end
xlim([1,zRange]);ylim([0,2000]);title('file1 Alines');ax=gca;ax.FontSize=18;

datB7 = zeros(size(datB5)) ;
for ii = 1:Dim.nx
    datB7(:,ii)=datB5(:,ii).*exp(1j.*dphase1);
end
RR0 = ifft(datB7,[],1);
RRB = abs(RR0(1:zRange,:,:));
figure
for ii = 2:50:Dim.nx
    plot(RRB(:,ii),'LineWidth',2);hold on
end
xlim([1,zRange]);ylim([0,2000]);title('file2 Alines');ax=gca;ax.FontSize=18;

%% phase stability check
figure;plot((2:Dim.nx).*16/1000,angle(RR0(128,2:end)),'LineWidth',2);title('file1 phase stability, std=0.054');ylabel('rad');xlabel('time(ms)');ax=gca;ax.FontSize=18;
std(angle(RR0(128,2:end)))

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% process all files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
for ifile = 1:length(files)
    filename=files(ifile).name;
    filePath=[datapath,filename];
    Dim.nk = 1024;
    % read spectrum
    fid=fopen(filePath,'r','l');
    dat = fread(fid, Dim.nk*Dim.nxRpt*Dim.nx*Dim.nyRpt, 'uint16');
    dat=reshape(dat, [Dim.nk Dim.nxRpt*Dim.nx Dim.nyRpt]);
    % trim non-interference samples
    dat=dat(101:700,:);
    %
    dat2 = zeros(size(dat)) ;
    for ii = 1:Dim.nx
        tmp = double(dat(:,ii)) -smooth(dat(:,ii),11);% DC removal
        dat2(:,ii)=get_interp(dphase,linPhase,tmp,indice).*exp(1j.*dphase1);% interpolation and dispersion
    end
    RR0 = ifft(dat2,[],1);
    RR = abs(RR0(1:zRange,:,:));
    plot(RR(:,round(Dim.nx/2)),'LineWidth',2);hold on
end
fid = fopen('intpX.bin','wb');
fwrite(fid,single(dphase));
fclose(fid);
fid = fopen('intpXP.bin','wb');
fwrite(fid,single(dphase));
fclose(fid);
fid = fopen('intpIndice.bin','wb');
fwrite(fid,single(dphase));
fclose(fid);
fid = fopen('dspPhase.bin','wb');
fwrite(fid,single(dphase));
fclose(fid);