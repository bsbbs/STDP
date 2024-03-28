%% setup
addpath(genpath('./utils/'));
addpath(genpath('./func/'));
pastelPalette = [ ...
    0.8941, 0.4471, 0.2471;  % Pastel Red
    0.3686, 0.6118, 0.8275;  % Pastel Blue
    0.4314, 0.8039, 0.6314;  % Pastel Green
    0.9569, 0.7608, 0.5176;  % Pastel Yellow
    0.8157, 0.5412, 0.7647;  % Pastel Purple
    1.0000, 0.6471, 0.5098;  % Pastel Orange
    ];
OKeeffe = [
    255, 192, 203;  % Pink (Pale Violet Red)
    100, 149, 237;  % Cornflower Blue
    255, 127, 80;   % Coral
    144, 238, 144;  % Light Green (Pale Green)
    255, 228, 196;  % Bisque
    147, 112, 219;  % Medium Purple
    0, 206, 209;    % Dark Turquoise
    250, 128, 114;  % Salmon
    152, 251, 152;  % Pale Green (Light Green)
    218, 112, 214;  % Orchid
    ]/255;
h = figure;
hold on;
for i = 1:size(OKeeffe,1)
    plot(i, 0, '.', 'Color', OKeeffe(i,:), 'MarkerSize',180);
end
mysavefig(h, 'Colorspace', gnrloutdir, 14, [8,3], 0);
close(h);