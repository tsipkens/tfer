
clear;
close all;
clc;

addpath autils;

% Will not work without cmap folder locally.
addpath cmap;
cm = turbo;


%%
%== DMA ==================================================================%
d = logspace(log10(5), log10(3e3), 500)';  % reconstruction points
d_star = logspace(log10(7), log10(800), 5);  % mobility setpoints
z = 0:3;

prop = prop_dma()  % default DMA properties
% prop.Q_a = prop.Q_a .* ones(size(d_star));  prop.Q_a(end-2:end) = 2.5e-5;  % test variable flow rate
% prop.Q_c = prop.Q_c .* ones(size(d_star));  prop.Q_c(end-2:end) = 2.5e-4;
% prop.Q_s = prop.Q_a;
% prop.Q_m = prop.Q_c

[Adma, ~, prop] = tfer_dma(d_star, d, z, prop);

% Triangular version of the transfer function.
% Only for singly charged.
Admat = tfer_tri(d_star, d, prop.Rd);

%-{
f1 = figure(1);
cmap_sweep(length(d_star), cm);
semilogx(d, sum(Adma, 3));
hold on;
semilogx(d, Admat, '--k');
hold off;
xlim([d(1), d(end)]); xlabel('d_m [nm]');
f1.Position(3) = 1000;
title('DMA');
%}


%%
%== PMA ==================================================================%
%   Currently assumes a one-to-one map m > dm. 
m = logspace(-3, 3, 2e3)';  % reconstruction points
m_star = logspace(-2.5, 2, 5)';  % mass-to-charge setpoints
z = 0:4;

prop = prop_pma;
prop = massmob.add(prop, 'soot')
d = (m .* 1e-18 ./ prop.m0) .^ (1 / prop.Dm);  % get mobility diameters

sp = get_setpoint(prop, 'm_star', m_star .* 1e-18, 'Rm', 3);
Af = tfer_pma(sp, m, d, z, prop, '1C_diff');  % unipolar/Fuchs

%-{
f2 = figure(2);
cmap_sweep(length(m_star), cm);
semilogx(m, sum(Af, 3));
xlim([m(1), m(end)]); xlabel('m_p [fg]');
f2.Position(2) = 600;
f2.Position(3) = 1000;
title('PMA');
%}

% Triangular version of the transfer function.
Aft = tfer_tri(sp, m .* 1e-18, prop.zet, z);

%-{
hold on;
semilogx(m, Aft(:,:,3), 'k--');
hold off;
%}

%%
%== CHARGING =============================================================%
d_star = logspace(log10(13.1), log10(500), 10)';  % mobility setpoints

z3 = -6:6;
[Ac3, ~, model] = charger(d_star, z3);  % unipolar/Fuchs

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
Ac4 = charger(d_star, z4, [], 'fuchs', opts);  % unipolar/Fuchs

%-{
subplot(5, 1, 2);
cmap_sweep(length(d_star), cm);
semilogx([5e-1,z4(2:end)], Ac4, 'o-');
xlabel(['Fuchs, nit = ', num2str(opts.nit/1e12), 'x10^{12}']);

opts.eps = 13.5;
opts.nit = 5e11;
z4 = 0:100;
Ac4 = charger(d_star, z4, [], 'fuchs', opts);  % unipolar/Fuchs

subplot(5, 1, 4);
cmap_sweep(length(d_star), cm);
semilogx([5e-1,z4(2:end)], Ac4, 'o-');
xlabel(['Fuchs, nit = ', num2str(opts.nit/1e12), 'x0^{12}']);
%}


opts.eps = 13.5;
opts.nit = 4e13;
z5 = 0:300;  % Li model currently requires contiguous charges starting at zero and must exceed 127
Ac5 = charger(d_star, z5, [], 'li', opts);  % unipolar/Fuchs

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
Ac5 = charger(d_star, z5, [], 'li', opts);  % unipolar/Fuchs

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


%%
%== AAC ==================================================================%
da = logspace(0.8, 3.2, 1.5e3)';  % reconstruction vector
da_star = logspace(1, 3, 10);
prop = prop_aac();
prop = massmob.add(prop, 'soot')

% Limited trajectrory variant.
opts = struct();  opts.model = 'lt';
Aa = tfer_aac(da_star, da, prop, opts);

% Particle streamline.
opts = struct();  opts.diffusion = false;
Aa2 = tfer_aac(da_star, da, prop, opts);

% Scanning version (limited trajectory only).
opts = struct();
opts.model = 'lt';
opts.scan = 1;
prop.tsc = 300;
prop.omega_s = 4e3;
prop.omega_e = 10;
Aas = tfer_aac(da_star, da, prop, opts);  % uses default elpi properties

%-{
f6 = figure(6);
cmap_sweep(length(da_star), cm);
semilogx(da, Aa);
hold on;
semilogx(da, Aa2, 'k--');  % particle streamline
semilogx(da, Aas, 'k:');
hold off;
xlim([da(1), da(end)]); xlabel('d_a [nm]');
f6.Position(2) = 400;
f6.Position(3) = 1000;
title('AAC');
%}

