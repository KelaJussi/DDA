close all
crnt = pwd();
cd("C:\Users\kelavuoj\OneDrive - TUNI.fi\Jussi\DDA\Data");

%{
load('Data3D_square_n_1_51Px_662Py_400_theta0_phi90_psi90_N50_50_1_tau0_2_lsprx500_lspry330_cs0_16.mat')
P = -Result.P(:,:,1);
%emt ei se ehkä ihan näin toimi. Mut 100x100 toimii
x = 100; y = 100;

Ptest = P(2:3:end);
figure
%plot(imag(Ptest(400*25+[1:50])))
hold on
%plot(real(Ptest(400*25+[1:50])))
plot(imag(Ptest(50:50:end)))
plot(real(Ptest(50:50:end)))
plot(abs(Ptest(25:50:end)))
legend(["imag","real","abs"])


load('Data3D_square_n_1_51Px_662Py_400_theta0_phi90_psi90_N400_50_1_tau0_2_lsprx500_lspry330_cs0_16.mat')
P = -Result.P(:,:,1);
%emt ei se ehkä ihan näin toimi. Mut 100x100 toimii


Ptest = P(2:3:end);
figure
%plot(imag(Ptest(400*25+[1:50])))
hold on
%plot(real(Ptest(400*25+[1:50])))
plot(imag(Ptest(25:50:end)))
%plot(real(Ptest(25:400:end)))
%plot(abs(Ptest(25:400:end)))
legend(["imag","real","abs"])
%}


x = 100; y = 100;
LSPRy = 330;
Px = 662;
 load(['Data3D_square_n_1_51Px_' num2str(Px) 'Py_400_theta1_phi0_psi0_N',num2str(x),'_',num2str(y),'_1_tau0_2_lsprx500_lspry' num2str(LSPRy) '_cs0_16.mat'])
P = Result.P(:,:,:);
P = permute(P, [2 1 3]);
wl = Result.Array.wl;
T = exp(-1e9*Result.C.ext);
%emt ei se ehkä ihan näin toimi. Mut 100x100 toimii
%Ei toimi, tuntuu että jos xxx!=yyy niin tulee tähän liittyvää kohinaa.
%esim 400x50 hila tuottaa 8 partsun periodein piikkejä. En ymmärrä...
%pitäisi olla (110,100) -> 100 ekaa partsut (0,0-100), 110 settiä näitä 
lambsize = size(P,3);
P = reshape(P,3,y*x,lambsize);
P = squeeze(sum(P(:,:,:),1));
%P = squeeze(vecnorm(P,2,1));

%P = fftshift(P,1);
%P = permute(P,[2 1 3]);
%
%ytemp = y; y = x; x = ytemp;

X = 1:x; Y = 1:y;
%
Pwave = zeros(x,y,lambsize);
for ii = Y
    Pwave(:,ii,:) = P([((ii-1)*x+1):1:ii*x],:);
end
%Pwave = fftshift(Pwave,2);
%Pwave = fftshift(Pwave,1);
%}
Ptotal = squeeze(vecnorm(P,2,1));
%Pwave = permute(reshape(P,[y,x,lambsize]),[1 2 3]);
%Pwave = P;
Pwave = pagetranspose(Pwave);
%figure
%plot(1:x, squeeze(abs(P(2,50,:))))
%figure
%plot(imag(Ptest(400*25+[1:400])))
%hold on
%plot(real(Ptest(400*25+[1:400])))
%plot(-imag(P1d(25:50:end)))
%plot(-real(P1d(25:50:end)))
%plot(abs(P1d(25:50:end)))
%legend(["imag","real","abs"])
figure
subplot(222)
Pwave2 = Pwave;
Pwave2(Pwave2>2e-22) = 2e-22; 
surf(wl*1e9,X,squeeze(abs(Pwave(:,ceil(y/2),:))),'EdgeColor','interp')
ylabel("\(x\) (a.u.)")
xlabel("Wavelength (nm)")
title("Dipolemoment magnitude, \(x\) (a.u.) to Wavelength (nm)")
title("\(||\mathbf{p}||\)")
colorbar
xlim([960,1040])
view(2)
%zlim(1e-22*[1.8,4])
yticks(0:25:100);
yticklabels(num2str([0:0.25:1].'))
subplot(223)
wl2 = 986;
plot(X,squeeze(abs(Pwave(:,ceil(y/2),ceil((1040-wl2)/0.4)))), 'linewidth', 2)
hold on 
wl2 = 1020.4;
plot(X,squeeze(abs(Pwave(:,ceil(y/2),ceil((1040-wl2)/0.4)))), 'linewidth', 2)

legend("\(\lambda = 986\) nm", "\(\lambda = 1020.4\) nm", "\(\lambda = 1015.6\) nm", "\(\lambda = 1016\) nm", "\(\lambda = 1021\) nm")
xlabel("\(x\) (a.u.) ")
ylabel("\(||\mathbf{p}||\)")
xticks(0:25:100);
xticklabels(num2str([0:0.25:1].'))
ylim([2e-22,5.2e-22])
box on 
%surf(wl*1e9,Y,squeeze(abs(Pwave(ceil(x/2),:,:))),'EdgeColor','interp')
%ylabel("y")
%xlabel("Wavelength")
%xticks(flip(wl(1:lambsize/5:lambsize))*1e9)
%title("Dipolemoments, y to wavelength")
subplot(224)
[~,index] = min(T);

surf(Y,X,squeeze(abs(Pwave(:,:,index))),'EdgeColor','interp')
view(2)
colorbar
yticks(0:25:100)
xticks(0:25:100)
xticklabels(num2str([0:0.25:1].'))
yticklabels(num2str([0:0.25:1].'))
vec = [0];
Pverho = squeeze(Pwave(:,ceil(y/2),index+vec));
Pverho = Pverho./max(Pverho);
save("DipoleMomentEnvelopeDDA", "Pverho")
%{
for ii = 1:length(vec)
    plot3(X,vec(ii)*ones(length(X)),Pverho(:,ii))
    hold on
end
%}
%legend(num2str([-2:2]))

ylabel("\(x\) (a.u.)")
xlabel(" \(y\) (a.u.)")
title("Dipolemoments in resonance wavelength")
title("\(||\mathbf{p}(\lambda=1020.4 \mathrm{nm})||\)")
%zlim(1e-22*[3,3.4])
subplot(221)
plot(wl*1e9,T, 'linewidth', 2)
hold on
%plot(wl*1e9,1-(max(T)-min(T))*(Ptotal./max(Ptotal)))
ylabel("Transmission")
xlabel("Wavelength (nm)")
%set(gcf,"Position",[200 200 1000 700])
sgtitle(['N_x = ', num2str(x), ', N_y = ', num2str(y), ', Pol= y'])
xlim([960,1040])
ylim([0.69,0.9])
%figure
%plot(,P(:,index))
%hold on 

%fontsize(gcf, 12,"points")

%heatmap(P(:,:,1))

cd(crnt);


