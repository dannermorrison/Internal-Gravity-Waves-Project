%% Parameters

clear
close all

% Number of discretized points in x-, z-direction
N_X = 128;
N_Y = 128;

% Time parameters
dt = 0.25; % time increment per iteration
T = 2000; % final time
T_draw = 25; % time interval at which plot is animated

% Viscosity
nu = 0.001;

%%% Spatial Variables and Wavenumbers

% 1. Real Space
dx = 2*pi/N_X;
dy = 2*pi/N_Y;
[x,y] = meshgrid((1:N_X)*dx, (1:N_Y)*dy);

% 2. Fourier Space
kx = [0:(N_X/2-1) (-N_X/2):-1];
ky = [0:(N_Y/2-1) (-N_Y/2):-1]';
[k, l] = meshgrid(kx, ky);

% Wavenumber Vector Norms
k_norm2 = k.^2 + l.^2; % norm-squared of wavenumber vector
k_norm2a = k_norm2; % non-zero norm-squared of wavenumber vector
k_norm2a(k_norm2a == 0) = 1; % prevents division by zero
%% Taylor-Green Vortex Example

% Initial Velocity
u = sin(x).*cos(y);
v = -cos(x).*sin(y);

%%% Simulation

for t = 0:dt:T

    % Vorticity
    zeta_k = 1i*k.*fft2(v) - 1i*l.*fft2(u); % Fourier Space

    % Dissipation
    zeta_k = zeta_k.*exp(-nu.*k_norm2.*dt);
    
    % Streamfunction
    psi_k = zeta_k./k_norm2a;
    u = real(ifft2(1i*l.*psi_k));
    v = real(ifft2(-1i*k.*psi_k)); 
    
    %%% Animation
    if (mod(t, T_draw) == 0)
        zeta = real(ifft2(zeta_k));
        pcolor(x,y,zeta); shading interp;
        axis square;
        colorbar; caxis([-2,2])
        hold on;
        quiver(x(1:8:N_X,1:8:N_Y),y(1:8:N_X,1:8:N_Y),...
            u(1:8:N_X,1:8:N_Y),v(1:8:N_X,1:8:N_Y),...
            0, 'r', 'LineWidth', 0.5, 'MaxHeadSize', 0.2); shading interp;
        hold off;
        drawnow
        
        % Plot Elements
        set(gca, 'FontName', 'Times', 'FontSize', 15)
        set([xlabel('x', 'Position', [3.2,-0.4]), ylabel('y')],...
            'FontName', 'Times', 'FontSize', 15, 'FontWeight', 'bold',...
            'Interpreter', 'latex')
    end    
        
end

%% Thesis Figures

% Initial Velocity
u = sin(x).*cos(y);
v = -cos(x).*sin(y);

save_times = 0:500:2000; % Times to produce plots

%%% Simulation

for t = 0:dt:T

    % Vorticity
    zeta_k = 1i*k.*fft2(v) - 1i*l.*fft2(u); % Fourier Space

    % Dissipation
    zeta_k = zeta_k.*exp(-nu.*k_norm2.*dt);
    
    % Streamfunction
    psi_k = zeta_k./k_norm2a;
    u = real(ifft2(1i*l.*psi_k));
    v = real(ifft2(-1i*k.*psi_k));
    
    % Save Image
    if ismember(t, save_times) == 1
        zeta = real(ifft2(zeta_k));
        pcolor(x,y,zeta); shading interp;
        axis square;
        colorbar; caxis([-2,2]);
        hold on;
        quiver(x(1:8:N_X,1:8:N_Y),y(1:8:N_X,1:8:N_Y),...
            u(1:8:N_X,1:8:N_Y),v(1:8:N_X,1:8:N_Y),...
            0, 'r', 'LineWidth', 0.8, 'MaxHeadSize', 0.8); shading interp;
        hold off;
        drawnow
        
        % Plot Elements
        set(gca, 'FontName', 'Times', 'FontSize', 15)
        set([xlabel('x', 'Position', [3.2,-0.4]), ylabel('y')], ...
            'FontName', 'Times', 'FontSize', 15, 'FontWeight', 'bold', ...
            'Interpreter', 'latex')
        
        % Export Plot
        filename = sprintf('time%d.png', t);
        exportgraphics(gcf, filename, 'Resolution', 400)
    end        
    
end

%% Additional Thesis Figures

% Initial Velocity
u = sin(x).*cos(y);
v = -cos(x).*sin(y);

% Streamfunction
psi_k = 1i*l.*fft2(u) - 1i*k.*fft2(v); % Fourier Space

save_times = 0:2000:2000; % Times to produce plots

%%% Simulation

for t = 0:dt:T

    % Vorticity
    zeta_k = 1i*k.*fft2(v) - 1i*l.*fft2(u); % Fourier Space

    % Dissipation
    zeta_k = zeta_k.*exp(-nu.*k_norm2.*dt);
    
    % Streamfunction
    psi_k = zeta_k./k_norm2a;
    u = real(ifft2(1i*l.*psi_k));
    v = real(ifft2(-1i*k.*psi_k));
    
    % Save Image
    if ismember(t, save_times) == 1
        zeta = real(ifft2(zeta_k));
        pcolor(x,y,zeta); shading interp;
        axis square;
        colorbar; caxis([-3,3])
        hold on;
        quiver(x(1:8:N_X,1:8:N_Y),y(1:8:N_X,1:8:N_Y),...
            u(1:8:N_X,1:8:N_Y),v(1:8:N_X,1:8:N_Y),...
            1, 'r', 'LineWidth', 0.8, 'MaxHeadSize', 0.8); shading interp;
        hold off;
        drawnow
        
        % Plot Elements
        set(gca, 'FontName', 'Times', 'FontSize', 15)
        set([xlabel('x', 'Position', [3.2,-0.4]), ylabel('y')], ...
            'FontName', 'Times', 'FontSize', 15, 'FontWeight', 'bold', ...
            'Interpreter', 'latex')
        
        % Export Plot
        filename = sprintf('AutoScale_time%d.png', t);
        exportgraphics(gcf, filename, 'Resolution', 400)
    end        
    
end