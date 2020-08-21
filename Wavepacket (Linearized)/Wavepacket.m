%% Parameters

clear
close all

% Number of discretized points in x-, z-direction
N_X = 128;
N_Z = 128;

% Scaling
N = 0.1; % Buoyancy frequency (time scale)
L = 1.5; % Horizontal length scale
H = 3; % Vertical length scale
delta = L/H;

% Time parameters
dt = 0.025; % time increment per iteration (scaled)
T = 200; % final time
T_draw = 25; % time interval at which plot is animated

% Viscosity
nu = 1e-10;

%%% Spatial Variables and Wavenumbers

% 1. Real Space (Nondimensionalized)
dx = 2*pi/N_X;
dz = 2*pi/N_Z;
[x,z] = meshgrid((1:N_X)*dx, (1:N_Z)*dz);

% 2. Fourier Space (Nondimensionalized)
kx = [0:(N_X/2-1) (-N_X/2):-1];
kz = [0:(N_Z/2-1) (-N_Z/2):-1]';
[k, m] = meshgrid(kx, kz);

% Wavenumber Vector Norms
k_norm2 = k.^2 + delta^2*m.^2; % norm-squared of wavenumber vector
k_norm2a = k_norm2; % non-zero norm-squared of wavenumber vector
k_norm2a(k_norm2a == 0) = 1; % prevents division by zero

%% Non-Dimensionalized Wavepacket Example
clf

% Parameters
B = 0.05; % wave amplitude
k0 = 25; % wavenumber x-direction
m0 = 25; % wavenumber z-direction
omega = sqrt(k0^2/(k0^2 + delta^2*m0^2)); % angular velocity
sigma = 0.6; % width of wavepacket envelope

% Horizontal vorticity
zeta_pack = @(t) (B*k0/omega)*cos(k0*x + m0*z - omega*t).*exp(-0.5*((x-pi).^2 + (z-pi).^2)/sigma^2); % exact solution (unscaled)
zeta_k = fft2(zeta_pack(0)); % Fourier space (nondimensionalized)

% Bouyancy acceleration
b_pack = @(t) B*cos(k0*x + m0*z - omega*t).*exp(-0.5*((x-pi).^2 + (z-pi).^2)/sigma^2); % exact solution (unscaled)
b_k = fft2(b_pack(0)); % Fourier space (nondimensionalized)

%%% Simulation

for t = 0:dt:T
    psi_k = zeta_k./k_norm2a;
    
    zeta_t{1} = -1i*k.*b_k; 
    b_t{1} = -1i*k.*psi_k;
    
    if t == 0 % Euler's Method
        zeta_k = zeta_k + dt*zeta_t{1};
        b_k = b_k + dt*b_t{1};
        
        zeta_t{2} = zeta_t{1};
        b_t{2} = b_t{1};       
        
    elseif t == dt % 2nd-Order Adams-Bashforth Method
        zeta_k = zeta_k + dt*0.5*(3*zeta_t{1} - zeta_t{2});
        b_k = b_k + dt*0.5*(3*b_t{1} - b_t{2});
        
        [zeta_t{2:3}] = zeta_t{1:2};
        [b_t{2:3}] = b_t{1:2};
        
    else % 3rd-Order Adams-Bashforth Method
        zeta_k = zeta_k + dt*0.083*(23*zeta_t{1} - 16*zeta_t{2} + 5*zeta_t{3});
        b_k = b_k + dt*0.083*(23*b_t{1} - 16*b_t{2} + 5*b_t{3});
        
        [zeta_t{2:3}] = zeta_t{1:2};
        [b_t{2:3}] = b_t{1:2};
        
    end
    
    % Dissipation
    zeta_k = zeta_k.*exp(-nu.*k_norm2.^2.*dt);
    b_k = b_k.*exp(-nu.*k_norm2.^2.*dt);
    
    %%% Animation
    if (mod(t, T_draw) == 0)
        zeta = real(ifft2(zeta_k));
        b = real(ifft2(b_k));
        subplot(1,2,1)
        pcolor(x,z,zeta); shading interp;
        axis square;
        colorbar;
        title('Numerical $\zeta(x,z,t)$', 'Interpreter', 'latex')
        subplot(1,2,2)
        pcolor(x,z,b); shading interp;
        axis square;
        colorbar;
        title('Numerical $b(x,z,t)$', 'Interpreter', 'latex')
        drawnow
    end    
end

%% Thesis Figures
clf

% Parameters
B = 0.05; % wave amplitude
k0 = 25; % wavenumber x-direction
m0 = 25; % wavenumber z-direction
omega = sqrt(k0^2/(k0^2 + delta^2*m0^2)); % angular velocity
sigma = 0.6; % width of wavepacket envelope

c_g = delta^2*m0*[m0, -k0]/(k0^2 + delta^2*m0^2)^1.5;

% Horizontal vorticity
zeta_pack = @(t) (B*k0/omega)*cos(k0*x + m0*z - omega*t).*exp(-0.5*((x-pi).^2 + (z-pi).^2)/sigma^2);
zeta_k = fft2(zeta_pack(0)); % Fourier space (nondimensionalized)

% Bouyancy acceleration
b_pack = @(t) B*cos(k0*x + m0*z - omega*t).*exp(-0.5*((x-pi).^2 + (z-pi).^2)/sigma^2);
b_k = fft2(b_pack(0)); % Fourier space (nondimensionalized)

save_times = [0,T]; % times to produce plots

%%% Simulation

for t = 0:dt:T
    psi_k = zeta_k./k_norm2a;
    
    zeta_t{1} = -1i*k.*b_k; 
    b_t{1} = -1i*k.*psi_k;
    
    if t == 0 % Euler's Method
        zeta_k = zeta_k + dt*zeta_t{1};
        b_k = b_k + dt*b_t{1};
        
        zeta_t{2} = zeta_t{1};
        b_t{2} = b_t{1};       
        
    elseif t == dt % 2nd-Order Adams-Bashforth Method
        zeta_k = zeta_k + dt*0.5*(3*zeta_t{1} - zeta_t{2});
        b_k = b_k + dt*0.5*(3*b_t{1} - b_t{2});
        
        [zeta_t{2:3}] = zeta_t{1:2};
        [b_t{2:3}] = b_t{1:2};
        
    else % 3rd-Order Adams-Bashforth Method
        zeta_k = zeta_k + dt*0.083*(23*zeta_t{1} - 16*zeta_t{2} + 5*zeta_t{3});
        b_k = b_k + dt*0.083*(23*b_t{1} - 16*b_t{2} + 5*b_t{3});
        
        [zeta_t{2:3}] = zeta_t{1:2};
        [b_t{2:3}] = b_t{1:2};
        
    end
    
    % Dissipation
    zeta_k = zeta_k.*exp(-nu.*k_norm2.^2.*dt);
    b_k = b_k.*exp(-nu.*k_norm2.^2.*dt);
    
    % Save Image
    if  t == save_times(1)
        
        % Vorticity
        zeta = real(ifft2(zeta_k));
        pcolor(x,z,zeta); shading interp;
        axis square;
        colorbar; caxis([-1.3,1.3]);
        hold on;
        scatter(x(1,64), z(64,1), 350, 'r.');
        hold off;
        drawnow;
        
        % Plot Elements
        set(gca, 'FontName', 'Times', 'FontSize', 15)
        set([xlabel('$x$*', 'Position', [3.2,-0.4]), ylabel('$z$*')], ...
            'FontName', 'Times', 'FontSize', 15, 'FontWeight', 'bold', ...
            'Interpreter', 'latex')
        
        % Export Plot
        filename = sprintf('zeta_time%d.png', t);
        exportgraphics(gcf, filename, 'Resolution', 400)
        
        % Buoyancy
        b = real(ifft2(b_k));
        pcolor(x,z,b); shading interp;
        axis square;
        colorbar; caxis([-0.04,0.04]);
        hold on;
        scatter(x(1,64), z(64,1), 350, 'r.');
        hold off;
        drawnow
        
        % Plot Elements
        set(gca, 'FontName', 'Times', 'FontSize', 15)
        set([xlabel('$x$*', 'Position', [3.2,-0.4]), ylabel('$z$*')], ...
            'FontName', 'Times', 'FontSize', 15, 'FontWeight', 'bold', ...
            'Interpreter', 'latex')
        
       % Export Plot
       filename = sprintf('b_time%d.png', t);
       exportgraphics(gcf, filename, 'Resolution', 400)

    elseif t == save_times(2)
        
        % Vorticity Plots
        
        zeta = real(ifft2(zeta_k));
        pcolor(x,z,zeta); shading interp;
        axis square;
        colorbar; caxis([-1.3,1.3]);
        hold on;
        scatter(x(1,64), z(64,1), 350, 'r.')
        quiver(x(1,64), z(64,1), T*c_g(1), T*c_g(2), 0,...
            'r', 'LineWidth', 1.5, 'MaxHeadSize', 0.5)
        scatter(x(1,64)+ T*c_g(1), z(64,1) + T*c_g(2), 350,'b.')
        hold off;
        drawnow;
        
        % Plot Elements
        set(gca, 'FontName', 'Times', 'FontSize', 15)
        set([xlabel('$x$*', 'Position', [3.2,-0.4]), ylabel('$z$*')], ...
            'FontName', 'Times', 'FontSize', 15, 'FontWeight', 'bold', ...
            'Interpreter', 'latex')
        
        % Export Plot
        filename = sprintf('zeta_time%d.png', t);
        exportgraphics(gcf, filename, 'Resolution', 400)
        
        % Buoyancy
        b = real(ifft2(b_k));
        pcolor(x,z,b); shading interp;
        axis square;
        colorbar; caxis([-0.04,0.04]);
        hold on;
        scatter(x(1,64), z(64,1), 350, 'r.')
        quiver(x(1,64), z(64,1), T*c_g(1), T*c_g(2), 0,...
            'r', 'LineWidth', 1.5, 'MaxHeadSize', 0.5)
        scatter(x(1,64)+ T*c_g(1), z(64,1) + T*c_g(2), 350,'b.')
        hold off;
        drawnow
        
        % Plot Elements
        set(gca, 'FontName', 'Times', 'FontSize', 15)
        set([xlabel('$x$*', 'Position', [3.2,-0.4]), ylabel('$z$*')], ...
            'FontName', 'Times', 'FontSize', 15, 'FontWeight', 'bold', ...
            'Interpreter', 'latex')
        
       % Export Plot
       filename = sprintf('b_time%d.png', t);
       exportgraphics(gcf, filename, 'Resolution', 400)
        
    end
end