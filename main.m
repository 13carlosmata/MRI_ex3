%% Simulation Exercise number3
addpath 'seemri'
clear all
close all
%%      1 Preparations
gammabar = 42.58e6;
gamma = 2*pi*gammabar;
%%      2 Slice Selection
B0 = 3;
iv1 = ImagingVolume(-1:0.25:5, -1:0.25:1, 0.8, 0.07, 1, 'PlotScale', 2);

tp = 2e-3;  B1 = 2.9e-6;
rf1 = SincPulse(B1, gammabar*B0, 0, tp);

g1 = Gradient([0 tp 1.5*tp], [70e-6 -70e-6 0]);
[S, ts] = seemri(iv1, B0, rf1, g1, [], ADC(1.5*tp, tp/100));
%%      2.1
rf2 = SincPulse(B1, gammabar*B0+4e3, 0, tp); rf3 = SincPulse(B1, gammabar*B0+8e3, 0, tp);
iv2 = ImagingVolume(-1:0.25:5, -1:0.25:1, 0.8, 0.07, 1, 'PlotScale', 2);
iv3 = ImagingVolume(-1:0.25:5, -1:0.25:1, 0.8, 0.07, 1, 'PlotScale', 2);
[S, ts] = seemri(iv2, B0, rf2, g1, [], ADC(1.5*tp, tp/100));
[S, ts] = seemri(iv3, B0, rf3, g1, [], ADC(1.5*tp, tp/100));
figure;
subplot(3,1,1); plot(iv1);
subplot(3,1,2); plot(iv2);
subplot(3,1,3); plot(iv3);
%%      2.2
g2 = Gradient([0 tp 1.5*tp], [110e-6 -110e-6 0]);
g3 = Gradient([0 tp 1.5*tp], [150e-6 -150e-6 0]);
iv1.toEquilibrium();
iv2.toEquilibrium();
iv3.toEquilibrium();
[S, ts] = seemri(iv1, B0, rf3, g1, [], ADC(1.5*tp, tp/100));
[S, ts] = seemri(iv2, B0, rf3, g2, [], ADC(1.5*tp, tp/100));
[S, ts] = seemri(iv3, B0, rf3, g3, [], ADC(1.5*tp, tp/100));
figure;
subplot(3,1,1); plot(iv1);
subplot(3,1,2); plot(iv2);
subplot(3,1,3); plot(iv3);
%%      2.3
tp = 2e-3; figure;
Delta_f = 12/tp;
Delta_z=2;
z0=0.001;
Gss=Delta_f/(gammabar*Delta_z);
g = Gradient([0 tp 1.5*tp], [Gss -Gss 0]);
f_rf=gammabar*(B0+Gss*Delta_z);
rf = SincPulse(B1, f_rf, 0, tp);
iv1.toEquilibrium();
[S, ts] = seemri(iv1, B0, rf, g, [], ADC(1.5*tp, tp/100));
%%
figure; iv1.toEquilibrium();
[S, ts] = seemri(iv1, B0, rf, g, g, ADC(1.5*tp, tp/100));
Gx=Gss;
Gss_N=(Gx^2+Gx^2)^0.5;
Delta_s=(Gx/Gss_N)*Delta_z;
%%
figure; u = @(x,y) sqrt(x.^2+y.^2)<=3;
x = -5:0.1:5; y = -5:0.1:5;
[xs, ys] = meshgrid(x, y); imagesc(x, y, u(xs, ys));
axis image; colormap gray; title('Image space')
Fu = @(kx,ky) 3*besselj(1, 2*pi*(sqrt(kx.^2+ky.^2)+1e-9*(kx==0 & ky==0))*3)...
./(sqrt(kx.^2+ky.^2)+1e-9*(kx==0 & ky==0));
kmax = 10; dk = 0.1;
ks = -kmax:dk:kmax-dk;
[kxs, kys] = meshgrid(ks,ks);
imagesc(ks, ks, Fu(kxs, kys));
axis image; colormap gray; title('k-space')
mrireconstruct(Fu(kxs, kys), kmax, 'Plot', true)
title(sprintf('kmax = %g, dk = %g', kmax, dk))

