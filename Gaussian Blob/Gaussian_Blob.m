%% Parameters

clear
close all

% Number of discretized points in x-, z-direction
N_X = 128;
N_Z = 128;

% Time parameters
dt = 0.025; % time increment per iteration (scaled)
T = 5; % final time
T_draw = 0.1; % time interval at which plot is animated

% Scaling
N = 0.1; % Buoyancy frequency (time scale)
L = 1.5; % Horizontal length scale
H = 3; % Vertical length scale
delta = L/H;

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

%% Gaussian Blob Example

% Time parameters
dt = 0.025; % time increment per iteration (scaled)
T = 5; % final time
T_draw = 0.1; % time interval at which plot is animated

% Viscosity
nu = 3e-8;

% Parameters
sigma = 0.4; % width of wavepacket envelope

% Horizontal vorticity
zeta_k = fft2(zeros(N_X, N_Z)); % Fourier space (nondimensionalized)

% Bouyancy acceleration
b_k = -fft2(2*ones(N_X,N_Z).*exp(-0.5*((x-pi).^2 + (z-pi).^2)/sigma^2)); % Fourier space (nondimensionalized)

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
        zeta = real(ifft2(zeta_k));
        b = real(ifft2(b_k));
%         subplot(1,2,1)
%         pcolor(x,z,zeta); shading interp;
%         axis square;
%         colorbar; caxis([-5,5]);
%         hold on;
%         quiver(x(1:8:N_X,1:8:N_Z),z(1:8:N_X,1:8:N_Z),...
%             u(1:8:N_X,1:8:N_Z),w(1:8:N_X,1:8:N_Z),...
%             0, 'r', 'LineWidth', 0.5, 'MaxHeadSize', 0.2); shading interp;
%         hold off;
%         title('Numerical $\zeta(x,z,t)$', 'Interpreter', 'latex')
%         subplot(1,2,2)
%         pcolor(x,z,b); shading interp;
%         axis square;
%         colorbar; caxis([-2,0]);
%                 hold on;
%         quiver(x(1:8:N_X,1:8:N_Z),z(1:8:N_X,1:8:N_Z),...
%             u(1:8:N_X,1:8:N_Z),w(1:8:N_X,1:8:N_Z),...
%             0, 'r', 'LineWidth', 0.5, 'MaxHeadSize', 0.2); shading interp;
%         hold off;
%         title('Numerical $b(x,z,t)$', 'Interpreter', 'latex')
%         drawnow
        
    end    
end
%% Thesis Figures #1 (Colorbar Scaling)
clf

% Viscosity
nu = 3e-8;

% Parameters
sigma = 0.4; % width of wavepacket envelope

save_times = 0:0.5:T; % Times to produce plots

% Horizontal vorticity
zeta_k = fft2(zeros(N_X, N_Z)); % Fourier space (nondimensionalized)

% Bouyancy acceleration
b_k = -fft2(2*ones(N_X,N_Z).*exp(-0.5*((x-pi).^2 + (z-pi).^2)/sigma^2)); % Fourier space (nondimensionalized)

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
    if ismember(t, save_times) == 1
        
        % Vorticity
        zeta = real(ifft2(zeta_k));
        pcolor(x,z,zeta); shading interp;
        axis square;
        colorbar; caxis([-5,5]);
        hold on;
        quiver(x(1:8:N_X,1:8:N_Z),z(1:8:N_X,1:8:N_Z),...
            u(1:8:N_X,1:8:N_Z),w(1:8:N_X,1:8:N_Z),...
            0, 'r', 'LineWidth', 0.5, 'MaxHeadSize', 0.2); shading interp;
        hold off;

        % Plot Elements
        set(gca, 'FontName', 'Times', 'FontSize', 15)
        set([xlabel('$x$', 'Position', [3.2,-0.4]), ylabel('$z$')], ...
            'FontName', 'Times', 'FontSize', 15, 'FontWeight', 'bold', ...
            'Interpreter', 'latex')
        
        % Export Plot
        filename = sprintf('zeta_time_scaled%d.png', t);
        exportgraphics(gcf, filename, 'Resolution', 400)
        
        % Buoyancy
        b = real(ifft2(b_k));
        pcolor(x,z,b); shading interp;
        axis square;
        colorbar; caxis([-2,2]);
        hold on;
        quiver(x(1:8:N_X,1:8:N_Z),z(1:8:N_X,1:8:N_Z),...
            u(1:8:N_X,1:8:N_Z),w(1:8:N_X,1:8:N_Z),...
            0, 'r', 'LineWidth', 0.5, 'MaxHeadSize', 0.2); shading interp;
        hold off;
        drawnow;
        
        % Plot Elements
        set(gca, 'FontName', 'Times', 'FontSize', 15)
        set([xlabel('$x$', 'Position', [3.2,-0.4]), ylabel('z')],...
            'FontName', 'Times', 'FontSize', 15, 'FontWeight', 'bold',...
            'Interpreter', 'latex')
        
        % Export Plot
        filename = sprintf('b_time_scaled%d.png', t);
        exportgraphics(gcf, filename, 'Resolution', 400)
    end    
end

%% Thesis Figures #2 (No Colorbar Scaling)
clf

% Viscosity
nu = 3e-8;

% Parameters
sigma = 0.4; % width of wavepacket envelope

save_times = 0:0.5:T; % Times to produce plots

% Horizontal vorticity
zeta_k = fft2(zeros(N_X, N_Z)); % Fourier space (nondimensionalized)

% Bouyancy acceleration
b_k = -fft2(2*ones(N_X,N_Z).*exp(-0.5*((x-pi).^2 + (z-pi).^2)/sigma^2)); % Fourier space (nondimensionalized)

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
    if ismember(t, save_times) == 1
        
        % Vorticity
        zeta = real(ifft2(zeta_k));
        pcolor(x,z,zeta); shading interp;
        axis square;
        colorbar;
        hold on;
        quiver(x(1:8:N_X,1:8:N_Z),z(1:8:N_X,1:8:N_Z),...
            u(1:8:N_X,1:8:N_Z),w(1:8:N_X,1:8:N_Z),...
            0, 'r', 'LineWidth', 0.5, 'MaxHeadSize', 0.2); shading interp;
        hold off;

        % Plot Elements
        set(gca, 'FontName', 'Times', 'FontSize', 15)
        set([xlabel('$x$', 'Position', [3.2,-0.4]), ylabel('$z$')], ...
            'FontName', 'Times', 'FontSize', 15, 'FontWeight', 'bold', ...
            'Interpreter', 'latex')
        
        % Export Plot
        filename = sprintf('zeta_time%d.png', t);
        exportgraphics(gcf, filename, 'Resolution', 400)
        
        % Buoyancy
        b = real(ifft2(b_k));
        pcolor(x,z,b); shading interp;
        axis square;
        colorbar;
        hold on;
        quiver(x(1:8:N_X,1:8:N_Z),z(1:8:N_X,1:8:N_Z),...
            u(1:8:N_X,1:8:N_Z),w(1:8:N_X,1:8:N_Z),...
            0, 'r', 'LineWidth', 0.5, 'MaxHeadSize', 0.2); shading interp;
        hold off;
        drawnow;
        
        % Plot Elements
        set(gca, 'FontName', 'Times', 'FontSize', 15)
        set([xlabel('$x$', 'Position', [3.2,-0.4]), ylabel('$z$')],...
            'FontName', 'Times', 'FontSize', 15, 'FontWeight', 'bold',...
            'Interpreter', 'latex')
        
        % Export Plot
        filename = sprintf('b_time%d.png', t);
        exportgraphics(gcf, filename, 'Resolution', 400)
    end    
end
