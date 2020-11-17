function[gassian_phi] = gaussDistribution(mode,SD,x_phi)
%% Gaussian distribution


x_mm = 2.^(-1.*x_phi);

mode_phi = mode;

a = (1/(SD*(sqrt(2*pi))));

gassian_phi = a.*exp(-0.5.*((x_phi - mode_phi)./SD).^2);


end
