clear 
close all

% Figure folder
figFolder = "figures/";

% Vector of tau values
Tau = 0:0.01:10;


% NDFE
B0 = zeros(size(Tau));

% Calculate the two BDFEs (for Tau > 4)
Bplus = nan(size(Tau));
Bminus = nan(size(Tau));
Bplus(Tau >= 4) = 0.5 * (1 + sqrt(1-4./Tau(Tau >= 4)));
Bminus(Tau >= 4) = 0.5 * (1 - sqrt(1-4./Tau(Tau >= 4)));

cols = colororder;

h = figure(1);
h.Position = [   100   502   738   392];
hold on
plot(Tau, B0, 'LineWidth', 2, 'Color', cols(1, :))
plot(Tau, Bplus, 'LineWidth', 2, 'Color', cols(2, :))
plot(Tau, Bminus, '--', 'LineWidth', 2, 'Color', cols(2, :))
xlabel('\tau')
ylabel('B^*')
ylim([0 1])
legend('NDFE', 'BDFE_+', 'BDFE_-', 'location', 'northwest')

ha = annotation('arrow');
ha.Parent = gca;
ha.X = 2 * [1 1]; 
ha.Y = [0.6, 0.4];
ha.Color = [0 0 0];

ha = annotation('arrow');
ha.Parent = gca;
ha.X = 7 * [1 1]; 
ha.Y = [0.15, 0.05];
ha.Color = [0 0 0];

ha = annotation('arrow');
ha.Parent = gca;
ha.X = 7 * [1 1]; 
ha.Y = [0.95, 0.85];
ha.Color = [0 0 0];

ha = annotation('arrow');
ha.Parent = gca;
ha.X = 7 * [1 1]; 
ha.Y = [0.45, 0.55];
ha.Color = [0 0 0];

fName = figFolder + "disease_free_bifurcation.png";
saveas(h, fName);

