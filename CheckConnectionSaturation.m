% N Exct connections every I received
tmp = sum(Ntwk.Cnnct_EI,2);
h = figure;
histogram(tmp);

% N Inhbt connections every E received
tmp = sum(Ntwk.Cnnct_IE,2);
h = figure;
histogram(tmp);

% N Exct connections every E received
tmp1 = sum(Ntwk.Cnnct_Input,2);
tmp2 = sum(Ntwk.Cnnct_EE,2);
h = figure;
histogram(tmp1+tmp2);

%% Summed Exct weights every E received from Inputs
tmp = sum(Ntwk.wInput,2);
h = figure;
histogram(tmp);