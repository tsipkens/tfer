
clear;
close all;
clc;

% addpath autils;

% Will not work without cmap folder locally.
addpath cmap;
cm = turbo;


%%
%== DMA ==================================================================%
d = logspace(log10(10), log10(800), 500)';  % reconstruction points
d_star = logspace(log10(13.1), log10(600), 20)';  % mobility setpoints
z = 0:4;

prop = prop_dma  % default DMA properties

Adma = tfer_dma(d_star' .* 1e-9, d .* 1e-9, z, prop);

% Triangular version of the transfer function.
Admat = tfer_tri(d_star', d, prop.Rd);

%-{
f1 = figure(1);
cmap_sweep(length(d_star), cm);
semilogx(d, Adma(:,:,2));
hold on;
semilogx(d, Admat, '--');
hold off;
xlim([d(1), d(end)]); xlabel('d_m [nm]');
f1.Position(3) = 1000;
title('DMA');
%}


%%
%== PMA ==================================================================%
%   Currently assumes a one-to-one map m > dm. 
m = logspace(-3, 2, 500)';  % reconstruction points
m_star = logspace(-2.5, 1.5, 20)';  % mass-to-charge setpoints
z = 0:100;  idx = 3;

prop = prop_pma;
prop = massmob.add(prop, 'soot')
d = (m .* 1e-18 ./ prop.m0) .^ (1 / prop.Dm);  % get mobility diameters

sp = get_setpoint(prop, 'm_star', m_star .* 1e-18, 'Rm', 3);
Af = tfer_pma(sp, m, d, z, prop, '1C_diff');  % unipolar/Fuchs

%-{
f2 = figure(2);
cmap_sweep(length(m_star), cm);
semilogx(m, Af(:,:,idx));
xlim([m(1), m(end)]); xlabel('m_p [fg]');
f2.Position(2) = 600;
f2.Position(3) = 1000;
title('PMA');
%}

% Triangular version of the transfer function.
Aft = tfer_tri(sp, m .* 1e-18, prop.zet, z);

%-{
hold on;
semilogx(m, Aft(:,:,idx), '--');
hold off;
%}

%%
%== CHARGING =============================================================%
d_star = logspace(log10(13.1), log10(500), 10)';  % mobility setpoints

z3 = -6:6;
[Ac3, ~, model] = charger(d_star .* 1e-9, z3);  % unipolar/Fuchs

%-{
f3 = figure(3);
f3.Position(2) = 200;
f3.Position(3) = 500;
f3.Position(4) = 600;

subplot(5, 1, 1);
cmap_sweep(length(d_star), cm);
plot(z3, Ac3, 'o-');
xlabel('Gopal.-Wieds.');
title('Charge distributions');
%}


opts.eps = 13.5;
opts.nit = 4e13;
z4 = 0:100;
Ac4 = charger(d_star .* 1e-9, z4, [], 'fuchs', opts);  % unipolar/Fuchs

%-{
subplot(5, 1, 2);
cmap_sweep(length(d_star), cm);
semilogx([5e-1,z4(2:end)], Ac4, 'o-');
xlabel(['Fuchs, nit = ', num2str(opts.nit/1e12), 'x10^{12}']);

opts.eps = 13.5;
opts.nit = 5e11;
z4 = 0:100;
Ac4 = charger(d_star .* 1e-9, z4, [], 'fuchs', opts);  % unipolar/Fuchs

subplot(5, 1, 4);
cmap_sweep(length(d_star), cm);
semilogx([5e-1,z4(2:end)], Ac4, 'o-');
xlabel(['Fuchs, nit = ', num2str(opts.nit/1e12), 'x0^{12}']);
%}


opts.eps = 13.5;
opts.nit = 4e13;
z5 = 0:300;  % Li model currently requires contiguous charges starting at zero and must exceed 127
Ac5 = charger(d_star .* 1e-9, z5, [], 'li', opts);  % unipolar/Fuchs

%-{
subplot(5, 1, 3);
cmap_sweep(length(d_star), cm);
semilogx([5e-1,z5(2:end)], Ac5, 'o-');
xlabel('z [-]');
xlim([5e-1, 100]);
xlabel(['Li, nit = ', num2str(opts.nit/1e12), 'x10^{12}']);

opts.eps = 13.5;
opts.nit = 5e11;
z5 = 0:300;  % Li model currently requires contiguous charges starting at zero and must exceed 127
Ac5 = charger(d_star .* 1e-9, z5, [], 'li', opts);  % unipolar/Fuchs

subplot(5, 1, 5);
cmap_sweep(length(d_star), cm);
semilogx([5e-1,z5(2:end)], Ac5, 'o-');
xlabel('z [-]');
xlim([5e-1, 100]);
xlabel(['Li, nit = ', num2str(opts.nit/1e12), 'x10^{12} / z [-]']);
%}


%%
%== BINNED DATA ==========================================================%
%   e.g., traditional SP2 interpretation. 
m = logspace(-3, 2, 500)';  % reconstruction vector
m_star = logspace(-2.5, 1.5, 20)';  % data vector

Abin = tfer_bin(m_star, m);

%-{
f4 = figure(4);
cmap_sweep(length(m_star), cm);
semilogx(m, Abin);
xlim([m(1), m(end)]);
f4.Position(2) = 500;
f4.Position(3) = 1000;
title('Binned');
%}


%%
%== ELPI/Impactor ========================================================%
da = logspace(0.5, 4.5, 500)';  % reconstruction vector

Aimp = tfer_elpi(da);  % uses default elpi properties

%-{
f5 = figure(5);
cmap_sweep(size(Aimp, 1), cm);
semilogx(da, Aimp);
xlim([da(1), da(end)]); xlabel('d_a [nm]');
f5.Position(2) = 400;
f5.Position(3) = 1000;
title('Impactor');
%}
