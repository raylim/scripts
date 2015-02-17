function plotNMF( prefix, nmfDir, minNumSig, maxNumSig )
% run NMF 
addpath(strcat(nmfDir, '/source/'));
addpath(strcat(nmfDir, '/plotting/'));
mkdir('temp');

minNumSig = str2num(minNumSig);
maxNumSig = str2num(maxNumSig);

for totalSignatures = minNumSig : maxNumSig
    inputFile = strcat(prefix, '_ts', num2str(totalSignatures), '.mat');
    load(inputFile);
    plotSignaturesToFile(prefix, processes, input, allProcesses, idx, processStabAvg);
end

quit
end
    
