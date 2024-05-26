function []=launcherscript(launcherparamfilename)
%Script to run PSO-based Matched Filtering on a data segment
%Inputs: segment number, txt file with filenames, txt file with segment
%details

%Load Files
fname = launcherparamfilename;
str = fileread(fname);
launcherparams = jsondecode(str);
hdf5filenames = launcherparams.hdf5filenames;
segfile = launcherparams.segfile;
psdfile = launcherparams.outfilename;
sampFreq = launcherparams.sampFreq;
segnum = launcherparams.segnum;
threshold = launcherparams.threshold;


%Load Data vector from file(s)
%Get Segment Data
segments = textread(segfile, '%s', 'delimiter', '\n');
Seg = split(segments{segnum});
segstart = str2num(Seg{2});
segend = str2num(Seg{3});
filenums = str2num(Seg{4});

%Get File data and load segment
if length(filenums) == 1
    files = textread(hdf5filenames, '%s', 'delimiter', '\n');
    fileStr = split(files{filenums});
    filename = fileStr{1};
    GPSstart = str2num(fileStr{2});
    
    %Get segment indices
    startidx = sampFreq*(segstart- GPSstart);
    endidx = sampFreq*(segend- GPSstart);

    %Load data from file and create segment data vector
    filedata = h5read(filename, '/strain/Strain')';
    segment = filedata(startidx+1:endidx);
    psdnum = filenums;
else
    segment = [];
    psdnum = filenums(1);
    for i = 1:length(filenums)
        files = textread(hdf5filenames, '%s', 'delimiter', '\n');
        fileStr = split(files{filenums(i)});
        filename = fileStr{1};
        GPSstart = str2num(fileStr{2});
        GPSend = str2num(fileStr{3});

        if GPSend > segstart && GPSstart < segstart && segend > GPSend
            %Get segment indices
            startidx = sampFreq*(segstart- GPSstart);
            endidx = sampFreq*(GPSend -GPSstart);
            %Load data from file and create segment data vector
            filedata = h5read(filename, '/strain/Strain')';
            seg = filedata(startidx:endidx);
            segment = [segment, seg];
        end
        if GPSstart < segend && GPSend > segend && segstart < GPSstart
            %Get segment indices
            startidx = 1;
            endidx = sampFreq*(segend - GPSstart);
            %Load data from file and create segment data vector
            filedata = h5read(filename, '/strain/Strain')';
            seg = filedata(startidx:endidx-1);
            segment = [segment, seg];
        end
    end
end

Seglen = length(segment)/sampFreq;

%{
signalinjectionsegs = [7,9,13,16,17,18];

%Inject signals
fmin = 30;
fmax = 700;
frange = [fmin,fmax];
initial_phase = 0; phase = 0;
siglen = Seglen;
if sum(ismember(signalinjectionsegs,segnum)) == 1
    %Draw parameters from uniform distribution
    m1 = unifrnd(1.4,15);
    m2 = unifrnd(1.4,15);
    masses = [m1,m2];
    %Estimate the corresponding chirptime parameters
    m1_val = m1*Msolar;
    m2_val = m2*Msolar;
    M = m1_val + m2_val;
    u = m1_val*m2_val/(m1_val + m2_val);
    n = u/M;
    tau0 = (5/(256*pi))*(1/fmin)*((G*M*pi*fmin/c^3)^(-5/3))*(1/n);
    tau1p5 = (1/8)*(1/fmin)*((G*M*pi*fmin/c^3)^(-2/3))*(1/n);

    r = unifrnd(30,100);
    mfac = unifrnd(2,3.2);
    ta = unifrnd(100,200);

    %Create signal
    signal = createsignal(siglen, frange, sampFreq, masses, r, initial_phase, phase, ta, mfac);

    %Record simulated strain signal in file
    sigparams = [segnum, tau0, tau1p5, m1, m2, ta, r, mfac];
    fileID = fopen('injectedsigs.txt','a');
    fprintf(fileID,'%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n',sigparams);
    fclose(fileID);
    %Add signal to segment
    segment = segment + signal;
end
%}

%Load Training Segment PSD
S = load(psdfile);
segPSDs = S.segPSDs;
trainidxs = S.trainidxs;
PSD = segPSDs{segnum};
trainidxs = trainidxs{segnum};


%Condition Data
rolloff = 0.125;
fmin = 30;
% paddedsegment = [zeros(1,sampFreq*rolloff), segment, zeros(1,sampFreq*rolloff)];
%Window and Highpass filter
segwin = segment.*tukeywin(length(segment),rolloff*sampFreq/length(segment))';
seghpass = highpass(segwin, fmin, sampFreq, ImpulseResponse="iir",Steepness=0.95);
filtsegment = seghpass;


%Whiten Segment
timeVec = (0:(length(filtsegment) - 1))/sampFreq;
[whtndseg,TFtotal] = segdatacond(filtsegment, PSD, sampFreq);


% Add pre-created, whitened custom signals in specific segments
%S = load('GVSsiginjections_final.mat'); sigsegnum = S.sigsegnums; signals = S.signals;
% S = load('GVSsiginjections_highsnr.mat'); sigsegnum = S.siginjsegs; signals = S.signals;
% index = find(sigsegnum == segnum);
% 
% if index
%     signal = signals(index,:);
%     whtndseg = whtndseg + signal;
% end
% 
% %Special modification for negative chirp-time space run on segment 878 which contains LIGO event GW190707_093326: Replaces the glitch near 347.8 seconds with WGN
% 
% if segnum == 878
%      whtndseg(347*sampFreq:348*sampFreq) = randn(1, length(whtndseg(347*sampFreq:348*sampFreq)));
% end


%Run PSO based matched-filtering on data vector
paramfile = 'allparamfiles.json';
runpso(segment, whtndseg, TFtotal, paramfile,segnum,Seglen);
end
