%% Examine if network initialization match the one trained.
% Conclusion is no problem
% Since the random generator are assigned with fixed number of random
% seeds, the networks generated will be the same based on the same code 
% RE the issue if network trained on HPC matches the local network
% structure
trainedfile = 'C:\Users\Bo\NYU Langone Health Dropbox\Shen Bo\Bo Shen Working files\STDP_Project\Analyses\ConstantValues\Sync_Val1.0\Weights_Seg14.mat';
initfile = 'C:\Users\Bo\NYU Langone Health Dropbox\Shen Bo\Bo Shen Working files\STDP_Project\Analyses\General1D\NtwkNMDA.mat';
load(trainedfile);
load(initfile);
% load another randomly generated network
tmp = load('C:\Users\Bo\NYU Langone Health Dropbox\Shen Bo\Bo Shen Working files\STDP_Project\Analyses\General1D\NtwkNMDA_II.mat');
NtwkII = tmp.Ntwk;

%% Compare the weights between two networks
h = figure;
imagesc(Ntwk.wInput - NtwkII.wInput);
sum(sum(Ntwk.wInput - NtwkII.wInput))
h = figure;
imagesc(Ntwk.wEI_initial - NtwkII.wEI_initial);
sum(sum(Ntwk.wEI_initial - NtwkII.wEI_initial))
h = figure;
imagesc(Ntwk.wEE_initial - NtwkII.wEE_initial);
sum(sum(Ntwk.wEE_initial - NtwkII.wEE_initial))
h = figure;
imagesc(Ntwk.Cnnct_Input - NtwkII.Cnnct_Input);
sum(sum(Ntwk.Cnnct_Input - NtwkII.Cnnct_Input))

%% Under the target network
visualize = 0;
NetworkgeneratorGPU1D;

EI = sum(WEI,1);
InputE = sum(Ntwk.wInput, 2);
h = figure;
subplot(2,2,1);
bar(InputE);
ylabel('Inputs received by E');
xlabel('E #');
subplot(2,2,2);
histogram(InputE);
xlabel('Inputs received by E');
ylabel('Frequency');
subplot(2,2,3);
bar(EI);
ylabel('Weights out from E');
xlabel('E #');
subplot(2,2,4);
histogram(EI)
xlabel('Weights out from E');
ylabel('Frequency');

h = figure;
subplot(2,1,1);
bar(EI(InputE < 2));
ylabel('E to I');
xlabel('Selected E(InputE < 2)');
subplot(2,1,2);
bar(InputE(EI < 3));
ylabel('Input to E');
xlabel('Selected E(EI < 3)');

CrossMtx = InputE*EI;
h = figure;
imagesc(CrossMtx);
