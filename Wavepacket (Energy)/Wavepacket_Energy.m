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

% De-aliasing Mask
M = ones(size(k_norm2));
M(abs(k) > (max(kx)*2/3)) = 0;
M(abs(m) > (max(kz)*2/3)) = 0;

%% Wavepacket Energy Example

% Parameters
B = 0.005; % wave amplitude
k0 = 25; % wavenumber x-direction
m0 = 25; % wavenumber z-direction
omega = sqrt(k0^2/(k0^2 + delta^2*m0^2)); % angular velocity
sigma = 0.6; % width of wavepacket envelope

% Horizontal vorticity
zeta_pack = @(t) (B*k0/omega)*cos(k0*x + m0*z - omega*t).*exp(-0.5*((x-pi).^2 + (z-pi).^2)/sigma^2);
zeta_k = fft2(zeta_pack(0)); % Fourier space (nondimensionalized)

% Bouyancy acceleration
b_pack = @(t) B*cos(k0*x + m0*z - omega*t).*exp(-0.5*((x-pi).^2 + (z-pi).^2)/sigma^2);
b_k = fft2(b_pack(0)); % Fourier space (nondimensionalized)

%%% Simulation

for t = 0:dt:T
    psi_k = zeta_k./k_norm2a;

    % Velocity
    u = real(ifft2(-1i*m.*psi_k));
    w = real(ifft2(1i*k.*psi_k));
    
    % Partials of zeta
    zeta_x = real(ifft2(1i*k.*zeta_k));
    zeta_z = real(ifft2(1i*m.*zeta_k));
    
    % Partials of b
    b_x = real(ifft2(1i*k.*b_k));
    b_z = real(ifft2(1i*m.*b_k));
    
    % Jacobian terms
    J_1 = w.*zeta_z + u.*zeta_x;
    J_2 = w.*b_z + u.*b_x;
    
    zeta_t{1} = -1i*k.*b_k - delta*M.*fft2(J_1); 
    b_t{1} = -1i*k.*psi_k - delta*M.*fft2(J_2);
    
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
        
        zeta = real(ifft2(zeta_k)); % vorticity
        b = real(ifft2(b_k)); % buoyancy
        
        % Energy
        E = 0.5*(u.^2 + delta^2*w.^2 + b.^2);
        [cx, cz] = find(E == max(max(E))); % coordinates of maximum energy

        pcolor(x,z,E); shading interp;
        axis square;
        colorbar; %caxis([0, 1.5e-3]);
        colormap jet;
        title('$E$($x$*, $z$*, $t$*)', 'Interpreter', 'Latex')
        xlabel('$x$*', 'fontweight','bold', 'Interpreter', 'Latex'); ylabel('$z_*$', 'fontweight','bold', 'Interpreter', 'Latex');
        hold on;
        scatter(x(1,cz), z(cx,1), 20,'wo', 'filled')
        hold off;
        drawnow;
    end 
end

%% Thesis Figures #1
clf

% Parameters
B = 0.05; % wave amplitude
k0 = 25; % wavenumber x-direction
m0 = 25; % wavenumber z-direction
omega = sqrt(k0^2/(k0^2 + delta^2*m0^2)); % angular velocity
sigma = 0.6; % width of wavepacket envelope

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

    % Velocity
    u = real(ifft2(-1i*m.*psi_k));
    w = real(ifft2(1i*k.*psi_k));
    
    % Partials of zeta
    zeta_x = real(ifft2(1i*k.*zeta_k));
    zeta_z = real(ifft2(1i*m.*zeta_k));
    
    % Partials of b
    b_x = real(ifft2(1i*k.*b_k));
    b_z = real(ifft2(1i*m.*b_k));
    
    % Jacobian terms
    J_1 = w.*zeta_z + u.*zeta_x;
    J_2 = w.*b_z + u.*b_x;
    
    zeta_t{1} = -1i*k.*b_k - delta*M.*fft2(J_1); 
    b_t{1} = -1i*k.*psi_k - delta*M.*fft2(J_2);
    
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
    if ismember(t, save_times) == 1
        
        % Vorticity
        zeta = real(ifft2(zeta_k));
        b = real(ifft2(b_k));

                % Energy
        E = 0.5*(u.^2 + delta^2*w.^2 + b.^2);
        [cx, cz] = find(E == max(max(E))); % coordinates of maximum energy

        pcolor(x,z,E); shading interp;
        axis square;
        colorbar; caxis([0,12e-4]);
        colormap jet;
        hold on;
        scatter(x(1,cz), z(cx,1), 30,'ko', 'filled')
        scatter(x(1,cz), z(cx,1), 20,'wo', 'filled')
        hold off;
        drawnow;
        
        % Plot Elements
        set(gca, 'FontName', 'Times', 'FontSize', 15)
        set([xlabel('$x$*', 'Position', [3.2,-0.4]), ylabel('$z$*')], ...
            'FontName', 'Times', 'FontSize', 15, 'FontWeight', 'bold', ...
            'Interpreter', 'latex')
        
        % Export Plot
        filename = sprintf('energy_time%d.png', t);
        exportgraphics(gcf, filename, 'Resolution', 400)
    end
end

%% Thesis Figures #2
clf

% Parameters
B = 0.005; % wave amplitude
k0 = 25; % wavenumber x-direction
m0 = 25; % wavenumber z-direction
omega = sqrt(k0^2/(k0^2 + delta^2*m0^2)); % angular velocity
sigma = 0.6; % width of wavepacket envelope

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

    % Velocity
    u = real(ifft2(-1i*m.*psi_k));
    w = real(ifft2(1i*k.*psi_k));
    
    % Partials of zeta
    zeta_x = real(ifft2(1i*k.*zeta_k));
    zeta_z = real(ifft2(1i*m.*zeta_k));
    
    % Partials of b
    b_x = real(ifft2(1i*k.*b_k));
    b_z = real(ifft2(1i*m.*b_k));
    
    % Jacobian terms
    J_1 = w.*zeta_z + u.*zeta_x;
    J_2 = w.*b_z + u.*b_x;
    
    zeta_t{1} = -1i*k.*b_k - delta*M.*fft2(J_1); 
    b_t{1} = -1i*k.*psi_k - delta*M.*fft2(J_2);
    
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
    if ismember(t, save_times) == 1
        
        % Vorticity
        zeta = real(ifft2(zeta_k));
        b = real(ifft2(b_k));

                % Energy
        E = 0.5*(u.^2 + delta^2*w.^2 + b.^2);
        [cx, cz] = find(E == max(max(E))); % coordinates of maximum energy

        pcolor(x,z,E); shading interp;
        axis square;
        colorbar; caxis([0,12e-6]);
        colormap jet;
        hold on;
        scatter(x(1,cz), z(cx,1), 30,'ko', 'filled')
        scatter(x(1,cz), z(cx,1), 20,'wo', 'filled')
        hold off;
        drawnow;
        
        % Plot Elements
        set(gca, 'FontName', 'Times', 'FontSize', 15)
        set([xlabel('$x$*', 'Position', [3.2,-0.4]), ylabel('$z$*')], ...
            'FontName', 'Times', 'FontSize', 15, 'FontWeight', 'bold', ...
            'Interpreter', 'latex')
        
        % Export Plot
        filename = sprintf('energy_time%da.png', t);
        exportgraphics(gcf, filename, 'Resolution', 400)
    end
end