%% Parameters

clear
close all

% Spatial Variables and Wavenumbers

% Number of discretized points in x-, z-direction
N_X = 128;
N_Z = 128;

% 1. Real Space (Nondimensionalized)
dx = 2*pi/N_X;
dz = 2*pi/N_Z;
[x,z] = meshgrid((1:N_X)*dx, (1:N_Z)*dz);

% 2. Fourier Space (Nondimensionalized)
kx = [0:(N_X/2-1) (-N_X/2):-1];
kz = [0:(N_Z/2-1) (-N_Z/2):-1]';
[k, m] = meshgrid(kx, kz);

% Wavenumber Vector Norms
k_norm2 = k.^2 + m.^2; % norm-squared of wavenumber vector
k_norm2a = k_norm2; % non-zero norm-squared of wavenumber vector
k_norm2a(k_norm2a == 0) = 1; % prevents division by zero

% De-aliasing Mask
M = ones(size(k_norm2));
M(abs(k) > (max(kx)*2/3)) = 0;
M(abs(m) > (max(kz)*2/3)) = 0;

%% Shear Flow Initialization #1
clf

% Horizontal Shear Flow
d = 0.5; U0 = 1;
U = U0*sech((z-pi)/d).^2;

N = 0.5; % Buoyancy frequency

% Richardson Number (Buoyancy Term/Flow Shear Term)
Ri = N^2*d^2*cosh((z-pi)/d).^4.*coth((z-pi)/d).^2/(4*U0^2);
disp('Richardson Number');
disp(min(Ri, [],'all'));
pcolor(x,z,U); shading interp;
axis square;
colorbar; cmap = parula(100);
colormap(cmap(50:100,:));
drawnow;

% Plot Elements
set(gca, 'FontName', 'Times', 'FontSize', 15)
set([xlabel('$x$', 'Position', [3.2,-0.4]), ylabel('$z$')], ...
    'FontName', 'Times', 'FontSize', 15, 'FontWeight', 'bold', ...
    'Interpreter', 'latex')

% Export Plot
filename = sprintf('shear_background1.png');
exportgraphics(gcf, filename, 'Resolution', 400)

%% Initial Conditions

colormap default;
% Horizontal vorticity
zeta_k = fft2(zeros(N_X, N_Z)); % Fourier space (nondimensionalized)

% Bouyancy acceleration
b_k = fft2(0.005*cos(2*x).*cos(z)); % Fourier space (nondimensionalized)

% Total Buoyancy
b = real(ifft2(b_k));
pcolor(x,z,N^2*z+b); shading interp;
axis square;
colorbar;
colormap = parula;

% Plot Elements
set(gca, 'FontName', 'Times', 'FontSize', 15)
set([xlabel('$x$', 'Position', [3.2,-0.4]), ylabel('$z$')], ...
    'FontName', 'Times', 'FontSize', 15, 'FontWeight', 'bold', ...
    'Interpreter', 'latex')

% Export Plot
filename = sprintf('totalb_init%d.png', 0);
exportgraphics(gcf, filename, 'Resolution', 400)

% Buoyancy
b = real(ifft2(b_k));
pcolor(x,z,b); shading interp;
axis square;
colorbar;
drawnow;

% Plot Elements
set(gca, 'FontName', 'Times', 'FontSize', 15)
set([xlabel('$x$', 'Position', [3.2,-0.4]), ylabel('$z$')], ...
    'FontName', 'Times', 'FontSize', 15, 'FontWeight', 'bold', ...
    'Interpreter', 'latex')

% Export Plot
filename = sprintf('b_init%d.png', 0);
exportgraphics(gcf, filename, 'Resolution', 400)


%% Shear Flow Example (Linearized)
clf;

% Time parameters
dt = 0.05; % time increment per iteration (scaled)
T = 7; % final time
T_draw = 0.1; % time interval at which plot is animated

% Viscosity
nu = 3e-4;

% Horizontal vorticity
zeta_k = fft2(zeros(N_X, N_Z)); % Fourier space (nondimensionalized)

% Bouyancy acceleration
b_k = fft2(0.005*cos(2*x).*cos(z)); % Fourier space (nondimensionalized)

% Horizontal background flow
U_k = fft2(U);

%%% Simulation

for t = 0:dt:T
    psi_k = zeta_k./k_norm2a;

    % Velocity
    u = real(ifft2(-1i*m.*psi_k));
    w = real(ifft2(1i*k.*psi_k));
    
    % Partials of U
    U_zz = real(ifft2(-m.^2.*U_k));
    
    % Partials of zeta
    zeta_x = real(ifft2(1i*k.*zeta_k));
    
    % Partials of b
    b_x = real(ifft2(1i*k.*b_k));
    
    % Nonlinear terms
    zeta_1 = fft2(U.*zeta_x);
    zeta_2 = fft2(U_zz.*w);
    b_1 = fft2(U.*b_x);
    
    zeta_t{1} = -1i*k.*b_k - M.*zeta_1 - M.*zeta_2;
    b_t{1} = -N^2*1i*k.*psi_k - M.*b_1;
    
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
        disp(t)

        % Vorticity
        zeta = real(ifft2(zeta_k));
        subplot(2,2,1)
        pcolor(x,z,zeta); shading interp;
        axis square;
        colorbar; %caxis([-0.01,0.01])
        title('$\zeta(x,z,t)$', 'Interpreter', 'Latex')
        xlabel('$x$', 'fontweight','bold', 'Interpreter', 'Latex'); ylabel('$z$', 'fontweight','bold', 'Interpreter', 'Latex')

        % Total Buoyancy
        b = real(ifft2(b_k));
        subplot(2,2,2)
        pcolor(x,z,N^2*z+b); shading interp;
        axis square;
        colorbar; %caxis([0,1.5]);
        title('$b(x,z,t)$', 'Interpreter', 'Latex')
        xlabel('$x$', 'fontweight','bold', 'Interpreter', 'Latex'); ylabel('$z$', 'fontweight','bold', 'Interpreter', 'Latex')
        
        % Velocity Vector Field
        subplot(2,2,3)
        quiver(x(1:8:N_X,1:8:N_Z),z(1:8:N_X,1:8:N_Z),u(1:8:N_X,1:8:N_Z),w(1:8:N_X,1:8:N_Z), 1, 'g', 'LineWidth', 1, 'MaxHeadSize', 1);
        axis square;
        title('Velocity Vector Field ($u$, $w$)', 'Interpreter', 'Latex')
        xlabel('$x$', 'fontweight','bold', 'Interpreter', 'Latex'); ylabel('$z$', 'fontweight','bold', 'Interpreter', 'Latex')
        
        % Buoyancy
        subplot(2,2,4)
        pcolor(x,z,b); shading interp;
        axis square;
        colorbar; %caxis([-5e-3,5e-3])
        title('$b(x,z,0)$', 'Interpreter', 'Latex')
        xlabel('$x$', 'fontweight','bold', 'Interpreter', 'Latex'); ylabel('$z$', 'fontweight','bold', 'Interpreter', 'Latex')        
        drawnow;
    end 
end

%% Thesis Figures #1
clf;

% Time parameters
dt = 0.05; % time increment per iteration (scaled)
T = 20; % final time
T_draw = 0.1; % time interval at which plot is animated

save_times = 1:1:20;

% Viscosity
nu = 3e-4;

% Horizontal vorticity
zeta_k = fft2(zeros(N_X, N_Z)); % Fourier space (nondimensionalized)

% Bouyancy acceleration
b_k = fft2(0.005*cos(2*x).*cos(z)); % Fourier space (nondimensionalized)

% Horizontal background flow
U_k = fft2(U);

%%% Simulation
for t = 0:dt:T
    psi_k = zeta_k./k_norm2a;

    % Velocity
    u = real(ifft2(-1i*m.*psi_k));
    w = real(ifft2(1i*k.*psi_k));
    
    % Partials of U
    U_zz = real(ifft2(-m.^2.*U_k));
    
    % Partials of zeta
    zeta_x = real(ifft2(1i*k.*zeta_k));
    
    % Partials of b
    b_x = real(ifft2(1i*k.*b_k));
    
    % Nonlinear terms
    zeta_1 = fft2(U.*zeta_x);
    zeta_2 = fft2(U_zz.*w);
    b_1 = fft2(U.*b_x);
    
    zeta_t{1} = -1i*k.*b_k - M.*zeta_1 - M.*zeta_2;
    b_t{1} = -N^2*1i*k.*psi_k - M.*b_1;
    
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
        colorbar; caxis([-0.01,0.01])
        
        % Plot Elements
        set(gca, 'FontName', 'Times', 'FontSize', 15)
        set([xlabel('$x$', 'Position', [3.2,-0.4]), ylabel('$z$')], ...
            'FontName', 'Times', 'FontSize', 15, 'FontWeight', 'bold', ...
            'Interpreter', 'latex')
        
        % Export Plot
        filename = sprintf('zeta_time%d.png', t);
        exportgraphics(gcf, filename, 'Resolution', 400)
        
        
        % Total Buoyancy
        b = real(ifft2(b_k));

        pcolor(x,z,N^2*z+b); shading interp;
        axis square;
        colorbar; caxis([0,1.5]);
        
        % Plot Elements
        set(gca, 'FontName', 'Times', 'FontSize', 15)
        set([xlabel('$x$', 'Position', [3.2,-0.4]), ylabel('$z$')], ...
            'FontName', 'Times', 'FontSize', 15, 'FontWeight', 'bold', ...
            'Interpreter', 'latex')
        
        % Export Plot
        filename = sprintf('totalb_time%d.png', t);
        exportgraphics(gcf, filename, 'Resolution', 400)

          
        % Buoyancy
        
        pcolor(x,z,b); shading interp;
        axis square;
        colorbar; caxis([-5e-3,5e-3])

        % Plot Elements
        set(gca, 'FontName', 'Times', 'FontSize', 15)
        set([xlabel('$x$', 'Position', [3.2,-0.4]), ylabel('$z$')], ...
            'FontName', 'Times', 'FontSize', 15, 'FontWeight', 'bold', ...
            'Interpreter', 'latex')
        
        % Export Plot
        filename = sprintf('b_time%d.png', t);
        exportgraphics(gcf, filename, 'Resolution', 400)

        
    end 
end
%% Thesis Figures #2
clf;

% Time parameters
dt = 0.05; % time increment per iteration (scaled)
T = 7; % final time
T_draw = 0.1; % time interval at which plot is animated

save_times = 7;

% Viscosity
nu = 3e-4;

% Horizontal vorticity
zeta_k = fft2(zeros(N_X, N_Z)); % Fourier space (nondimensionalized)

% Bouyancy acceleration
b_k = fft2(0.005*cos(2*x).*cos(z)); % Fourier space (nondimensionalized)

% Horizontal background flow
U_k = fft2(U);

%%% Simulation
for t = 0:dt:T
    psi_k = zeta_k./k_norm2a;

    % Velocity
    u = real(ifft2(-1i*m.*psi_k));
    w = real(ifft2(1i*k.*psi_k));
    
    % Partials of U
    U_zz = real(ifft2(-m.^2.*U_k));
    
    % Partials of zeta
    zeta_x = real(ifft2(1i*k.*zeta_k));
    
    % Partials of b
    b_x = real(ifft2(1i*k.*b_k));
    
    % Nonlinear terms
    zeta_1 = fft2(U.*zeta_x);
    zeta_2 = fft2(U_zz.*w);
    b_1 = fft2(U.*b_x);
    
    zeta_t{1} = -1i*k.*b_k - M.*zeta_1 - M.*zeta_2;
    b_t{1} = -N^2*1i*k.*psi_k - M.*b_1;
    
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
        quiver(x(1:8:N_X,1:8:N_Z),z(1:8:N_X,1:8:N_Z),u(1:8:N_X,1:8:N_Z),w(1:8:N_X,1:8:N_Z), 1, 'r', 'LineWidth', 1, 'MaxHeadSize', 1);
        hold off;
        
        % Plot Elements
        set(gca, 'FontName', 'Times', 'FontSize', 15)
        set([xlabel('$x$', 'Position', [3.2,-0.4]), ylabel('$z$')], ...
            'FontName', 'Times', 'FontSize', 15, 'FontWeight', 'bold', ...
            'Interpreter', 'latex')
        
        % Export Plot
        filename = sprintf('zeta_time%d_new.png', t);
        exportgraphics(gcf, filename, 'Resolution', 400)
        
        
        % Total Buoyancy
        b = real(ifft2(b_k));

        pcolor(x,z,N^2*z+b); shading interp;
        axis square;
        colorbar; caxis([0,1.5]);
        
        % Plot Elements
        set(gca, 'FontName', 'Times', 'FontSize', 15)
        set([xlabel('$x$', 'Position', [3.2,-0.4]), ylabel('$z$')], ...
            'FontName', 'Times', 'FontSize', 15, 'FontWeight', 'bold', ...
            'Interpreter', 'latex')
        
        % Export Plot
        filename = sprintf('totalb_time%d.png', t);
        exportgraphics(gcf, filename, 'Resolution', 400)
          
        % Buoyancy
        
        pcolor(x,z,b); shading interp;
        axis square;
        colorbar; caxis([-5e-3,5e-3])

        % Plot Elements
        set(gca, 'FontName', 'Times', 'FontSize', 15)
        set([xlabel('$x$', 'Position', [3.2,-0.4]), ylabel('$z$')], ...
            'FontName', 'Times', 'FontSize', 15, 'FontWeight', 'bold', ...
            'Interpreter', 'latex')
        
        % Export Plot
        filename = sprintf('b_time%d_new.png', t);
        exportgraphics(gcf, filename, 'Resolution', 400)

        
    end 
end
