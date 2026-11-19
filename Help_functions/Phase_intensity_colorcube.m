
phase = linspace(-1,1,200);
intensity = linspace(0,0.6,200);

[P,I] = meshgrid(phase,intensity);
%Vmin = 0.6;  % minimum brightness
%V = 1 - (1 - Vmin) * (mag_norm);
HSV = cat(3,...
    (P+pi)/(2*pi), ... % hue
    ones(size(P)), ... % saturation
    I);                % value
HSV = cat(3,...
    (P+1)/2, ... % hue
    I, ... % saturation
    1-I);   
RGB = hsv2rgb(HSV);
figure
imagesc(phase,intensity,RGB)
set(gca,'YDir','normal')
xlabel('Phase')
ylabel('Intensity')
