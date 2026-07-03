%Estructura HFB para caso de interacción USDB
clc;clear;
format long
m=24; % Dimension matrices, grados de libertad 
N0=8; %Nº de particulas

% ==============================================================================
% 1. CONSTRUCCIÓN DE LA BASE DESACOPLADA (m-scheme)
% ==============================================================================
states = [];
idx = 1;
% -1 para Protones, +1 para Neutrones
for two_mt = [-1, 1] 
    for orbit = [205, 1001, 203] % Orden de subcapas solicitado
        if orbit == 205,     two_j = 5;
        elseif orbit == 1001, two_j = 1;
        elseif orbit == 203,  two_j = 3;
        end
        
        for two_mj = -two_j:2:two_j
            states(idx).orbit = orbit;
            states(idx).two_j = two_j;
            states(idx).two_mj = two_mj;
            states(idx).two_mt = two_mt;
            idx = idx + 1;
        end
    end
end
m = length(states); % Debería ser exactamente 24

% ==============================================================================
% 2. LECTOR DEL FICHERO DE INTERACCIÓN (USDB)
% ==============================================================================
fid = fopen('usdb_baseacoplada.txt', 'r');
if fid == -1, error('No se pudo abrir el archivo de interacción.'); end

fgetl(fid); % Saltar cabecera de texto

% Leer configuración de órbitas e información general
orbit_info = fscanf(fid, '%d', 5);
num_orbits = orbit_info(2);
orbit_labels = orbit_info(3:5);

% Energías monoparticulares del archivo USDB (-3.92570, -3.2079, 2.1117)
spe_usdb = fscanf(fid, '%f', num_orbits); 
core_info = fscanf(fid, '%f', 4); % Info del core

% Inicializar contenedor para los elementos acoplados V(a,b,c,d,J,T)
V_coupled = zeros(3, 3, 3, 3, 6, 2); % Índices mapeados: 1->205, 2->1001, 3->203

orbit_map = @(x) find(x == [205, 1001, 203]);
orbit_two_j = @(x) (x==205)*5 + (x==1001)*1 + (x==203)*3;

while ~feof(fid)
    tok = fscanf(fid, '%d', 1);
    if isempty(tok), break; end
    if tok == 0
        peek1 = fscanf(fid, '%d', 1);
        if peek1 == 1
            a = fscanf(fid, '%d', 1); b = fscanf(fid, '%d', 1);
            c = fscanf(fid, '%d', 1); d = fscanf(fid, '%d', 1);
            j_min = fscanf(fid, '%d', 1); j_max = fscanf(fid, '%d', 1);
            
            n_j = j_max - j_min + 1;
            t0_vals = fscanf(fid, '%f', n_j);
            t1_vals = fscanf(fid, '%f', n_j);
            
            ia = orbit_map(a); ib = orbit_map(b);
            ic = orbit_map(c); id = orbit_map(d);
            
            ja2 = orbit_two_j(a); jb2 = orbit_two_j(b);
            jc2 = orbit_two_j(c); jd2 = orbit_two_j(d);
            
            for J = j_min:j_max
                v0 = t0_vals(J - j_min + 1);
                v1 = t1_vals(J - j_min + 1);
                
                for T = 0:1
                    if T == 0, val = v0; else, val = v1; end
                    
                    % Fases de permutación antisimétrica
                    phase_ba = (-1)^((ja2 + jb2)/2 - J +1 - T);
                    phase_dc = (-1)^((jc2 + jd2)/2 - J +1 - T);
                    
                    % Rellenar por simetrías y hermiticidad
                    V_coupled(ia, ib, ic, id, J+1, T+1) = val;
                    V_coupled(ib, ia, ic, id, J+1, T+1) = phase_ba * val;
                    V_coupled(ia, ib, id, ic, J+1, T+1) = phase_dc * val;
                    V_coupled(ib, ia, id, ic, J+1, T+1) = phase_ba * phase_dc * val;
                    
                    V_coupled(ic, id, ia, ib, J+1, T+1) = val;
                    V_coupled(ic, id, ib, ia, J+1, T+1) = phase_ba * val;
                    V_coupled(id, ic, ia, ib, J+1, T+1) = phase_dc * val;
                    V_coupled(id, ic, ib, ia, J+1, T+1) = phase_ba * phase_dc * val;
                end
            end
        end
    end
end
fclose(fid);

% ==============================================================================
% 3. PROYECCIÓN AL ESQUEMA-M (Construcción de v_barra)
% ==============================================================================
v_barra = zeros(m, m, m, m);

for alpha = 1:m
    st_a = states(alpha); ia = orbit_map(st_a.orbit);
    for beta = 1:m
        st_b = states(beta); ib = orbit_map(st_b.orbit);
        
        if alpha == beta, continue; end
        two_MJ = st_a.two_mj + st_b.two_mj;
        two_MT = st_a.two_mt + st_b.two_mt;
        
        % Factores de normalización de la base acoplada
        norm_ab = sqrt(1 + double(st_a.orbit == st_b.orbit));
        
        for gamma = 1:m
            st_c = states(gamma); ic = orbit_map(st_c.orbit);
            for delta = 1:m
                st_d = states(delta); id = orbit_map(st_d.orbit);
                
                if gamma == delta, continue; end
                if (st_c.two_mj + st_d.two_mj ~= two_MJ) || (st_c.two_mt + st_d.two_mt ~= two_MT)
                    continue;
                end
                
                norm_cd = sqrt(1 + double(st_c.orbit == st_d.orbit));
                
                sum_V = 0.0;
                j_min_allowed = abs(st_a.two_j - st_b.two_j)/2;
                j_max_allowed = (st_a.two_j + st_b.two_j)/2;
                
                for J = j_min_allowed:j_max_allowed
                    if abs(two_MJ) > 2*J, continue; end
                    for T = 0:1
                        if abs(two_MT) > 2*T, continue; end
                        
                        cg_J_ab = cg_doubled(st_a.two_j, st_a.two_mj, st_b.two_j, st_b.two_mj, 2*J, two_MJ);
                        if cg_J_ab == 0, continue; end
                        cg_J_cd = cg_doubled(st_c.two_j, st_c.two_mj, st_d.two_j, st_d.two_mj, 2*J, two_MJ);
                        if cg_J_cd == 0, continue; end
                        
                        cg_T_ab = cg_doubled(1, st_a.two_mt, 1, st_b.two_mt, 2*T, two_MT);
                        if cg_T_ab == 0, continue; end
                        cg_T_cd = cg_doubled(1, st_c.two_mt, 1, st_d.two_mt, 2*T, two_MT);
                        if cg_T_cd == 0, continue; end
                        
                        V_c = V_coupled(ia, ib, ic, id, J+1, T+1);
                        
                        % Factor combinatorio estándar de cambio de base antisimétrica
                        factor = norm_ab * norm_cd;
                        
                        sum_V = sum_V + factor * cg_J_ab * cg_T_ab * cg_J_cd * cg_T_cd * V_c;
                    end
                end
                v_barra(alpha, beta, gamma, delta) = sum_V;
            end
        end
    end
end


level = zeros(1, 24);
for i = 1:24
    if states(i).orbit == 205,     level(i) = spe_usdb(1);
    elseif states(i).orbit == 1001, level(i) = spe_usdb(2);
    elseif states(i).orbit == 203,  level(i) = spe_usdb(3);
    end
end

%Matriz compleja aleatoria
%rng(1); %Semilla
%Z0 = (randn(m) + 1i*randn(m));
Z0=randn(m);
%Imponer antisimetria
Z0 = Z0 - Z0.';

% construir U
M = eye(m) + Z0'*Z0;
U0 = inv(sqrtm(M));
V0 = Z0*U0;

% V0=[0 0 0 0 0 0;0 0 0 0 0 0;0 0 1 0 0 0;0 0 0 1 0 0;0 0 0 0 0 0;0 0 0 0 0 0]; %Particulas en el nivel 2
% U0=[1 0 0 0 0 0;0 1 0 0 0 0;0 0 0 0 0 0;0 0 0 0 0 0;0 0 0 0 1 0;0 0 0 0 0 1];


% construir W con estructura HFB
 W = [U0 conj(V0); V0 conj(U0)];
% 
% % comprobar unitariedad
 error = norm(W'*W - eye(2*m));
% 
% %Comprobar condiciones U y V
 cond1 = norm(U0'*U0 + V0'*V0)-1;
 cond2 = norm(U0.'*V0 + V0.'*U0);
% 
% %Construcción de las matriz de densidad (rho) y el tensor de paridad (kappa)
 

 rho=conj(V0)*V0.';
 kappa=conj(V0)*U0.';

% %Comprobación de sus propiedades

 prop1=rho'-rho; %Hermiticidad de la matriz de densidad (prop1==zeros(N))
 prop2=kappa.'+kappa; %Antisimetria del tensor de paridad

% Definición de Delta y Gamma
%Construccion gamma 

Gamma = zeros(m,m);
for a=1:m
for c=1:m
    
    sum_val = 0;
    
    for d=1:m
    for b=1:m
        
        sum_val = sum_val + v_barra(a,b,c,d)*rho(d,b);
        
    end
    end
    
    Gamma(a,c) = sum_val;
    
end
end

%Construccion Delta
Delta = zeros(m,m);
for a=1:m
for b=1:m
    
    sum_val = 0;
    
    for c=1:m
    for d=1:m
       
        sum_val = sum_val + v_barra(a,b,c,d)*kappa(c,d);
        
    end
    end
    
    Delta(a,b) = 0.5*sum_val;
    
end
end


% Energías monoparticulares 
eps = level;  

t = diag(level);
h = t + Gamma;

%Def de H20

H20 = U0' * h * conj(V0) ...
    + U0' * Delta * conj(U0) ...
    - V0' * h.' * conj(U0) ...
    - V0' * conj(Delta) * conj(V0);


%Def de N20

N20 = U0' * conj(V0) - V0' * conj(U0);

%Comprobacion antisimetria
check_H20 = norm(H20.' + H20); %¿=0?
check_N20 = norm(N20.' + N20); %¿=0?

%Cálculo de E_HFB

E1 = trace(t * rho);
E2 = 0.5 * trace(Gamma * rho);
E3 = -0.5 * trace(Delta * conj(kappa));
E_HFB = E1 + E2 + E3;
N=trace(rho); %Hay que ajustarlo al valor deseado
%ES REAL?
check_EReal =norm(imag(E_HFB)); %¿=0?


%% Ajuste puro del constraint (Número de partículas)
% --- Parámetros de control ---
max_iter = 10;    % Máximo de iteraciones
tol = 1e-6;        % Tolerancia de convergencia

fprintf('Iteración | Num. Partículas | Error N \n');
fprintf('------------------------------------------\n');

U = U0;
V = V0;

for k = 1:max_iter
    % 1. Recalcular las densidades actuales
    rho = conj(V) * V.';   
    N_calc=real(trace(rho));
    N20 = U' * conj(V) - V' * conj(U);

    % 2. Evaluar el error
    error_N = N0 - N_calc;

    % Criterio de parada temprano
    if abs(error_N) < tol
        fprintf('Convergencia alcanzada en la iteración %d con N=%d\n', k,N_calc);
        break;
    end
    % 3. Cálculo del paso (eta_N)  (Ec. (B.8))
    eta_N = error_N / real(trace(N20 * N20'));

    % 4. Dirección de actualización Z (Solo usamos N20)
    Z = eta_N * N20;
    Z = 0.5 * (Z - Z.'); % Asegurar que la matriz Z sea estrictamente antisimétrica 

    % 5. Reconstrucción de U y V (Preservando unitariedad mediante Cholesky)
    A0 = eye(m) + Z.' * conj(Z);
    L0 = chol(A0, 'lower');

    Uold = U;
    Vold = V;

    % Actualizacion U,V (ec. 25)
    U = (Uold + conj(Vold) * conj(Z)) * (inv(L0))';
    V = (Vold + conj(Uold) * conj(Z)) * (inv(L0))';

    % 6. Verificación de unitariedad Bogoliubov
    tol_unitarity = 1e-7;
    condic1 = norm(U'*U + V'*V - eye(m));
    condic2 = norm(U.'*V + V.'*U);

    if condic1 > tol_unitarity || condic2 > tol_unitarity
        warning('Iteración %d: ¡Pérdida de unitariedad detectada!', k);
        fprintf('  Error Condición 1: %.2e\n', condic1);
        fprintf('  Error Condición 2: %.2e\n', condic2);
    end
    % Mostrar progreso
    fprintf('%9d | %15.4f | %10.2e \n', k, N_calc, error_N);
end





%% Minimización de Energía por Gradiente con Paso Adaptativo y Constraint
% --- Parámetros de Control ---
max_iter_E = 1000;     % Iteraciones máximas del bucle de energía (externo)
max_iter_N = 10;      % Iteraciones máximas del bucle de partículas (interno)
eta_E = 0.05;         % ¡Paso inicial agresivo! Se auto-ajustará solo.
tol_E = 1e-6;         % Tolerancia para la convergencia en energía
tol_N = 1e-6;         % Tolerancia para el número de partículas


fprintf('Iter Ext | Energía HFB | Num. Partículas | Delta E  | Iter Int\n');
fprintf('---------------------------------------------------------------------------------------\n');

% --- Inicialización para la gráfica ---
historial_E = [];    % Almacenará los valores de energía
historial_iter = []; % Almacenará el número de iteración exitosa
E_actual=E_HFB;

for iter_E = 1:max_iter_E
    % =========================================================================
    % 1. CÁLCULO DE GRADIENTES EN EL PUNTO ACTUAL
    % =========================================================================
    rho = conj(V) * V.';
    kappa = conj(V)*U.';
    % Gamma actual

    Gamma = zeros(m,m);

    for a=1:m
        for c=1:m

            sum_val = 0;

            for d=1:m
                for b=1:m

                    sum_val = sum_val + v_barra(a,b,c,d)*rho(d,b);

                end
            end

            Gamma(a,c) = sum_val;

        end
    end

    % Delta actual
    Delta = zeros(m,m);

    for a=1:m

        for b=1:m

            sum_val = 0;

            for c=1:m
                for d=1:m

                    sum_val = sum_val + v_barra(a,b,c,d)*kappa(c,d);

                end
            end

            Delta(a,b) = 0.5*sum_val;

        end
    end
    

    N_final = real(trace(rho));
    h = t + Gamma;
    H20 = U' * h * conj(V) + U' * Delta * conj(U) - V' * h.' * conj(U) - V' * conj(Delta) * conj(V);
    N20 = U' * conj(V) - V' * conj(U); 
    E_old=E_actual;
    % =========================================================================
    % 2. PROYECCIÓN DEL MULTIPLICADOR DE LAGRANGE (LAMBDA)
    % =========================================================================
    denom_lambda = trace(N20 * N20');
    if denom_lambda > 1e-12
        lambda = trace(H20 * N20') / denom_lambda;
    else
        lambda = 0;
    end

    % =========================================================================
    % 3. PASO DE GRADIENTE PARA LA ENERGÍA
    % =========================================================================
    Z_E = -eta_E * (H20 - lambda * N20); %Direccion de gradiente
    Z_E = 0.5 * (Z_E - Z_E.'); %Unitariedad Z (forzada)

    A0 = eye(m) + Z_E.' * conj(Z_E);
    L0 = chol(A0, 'lower');

    U_old=U;
    V_old=V;

    U = (U_old + conj(V_old) * conj(Z_E)) * (inv(L0))';
    V = (V_old + conj(U_old) * conj(Z_E)) * (inv(L0))';

  
    % =========================================================================
    % 4. BUCLE INTERNO: RESTAURACIÓN ESTRICTA DEL NÚMERO DE PARTÍCULAS
    % =========================================================================
    sub_iter = 0; 
    for k_N = 1:max_iter_N
        rho_sub = conj(V) * V.';
        N_sub = real(trace(rho_sub));
        error_N = N0 - N_sub;

        if abs(error_N) < tol_N
            sub_iter = k_N - 1;
            break;
        end

        N20_sub = U' * conj(V) - V' * conj(U);
        norm_N20 = real(trace(N20_sub * N20_sub'));

        if norm_N20 < 1e-14, break; end

        eta_N = error_N / norm_N20;
        Z_N = eta_N * N20_sub;
        Z_N = 0.5 * (Z_N - Z_N.');

        A_sub = eye(m) + Z_N.' * conj(Z_N);
        L_sub = chol(A_sub, 'lower');

        Uold_sub = U; Vold_sub = V;
        U = (Uold_sub + conj(Vold_sub) * conj(Z_N)) * (inv(L_sub))';
        V = (Vold_sub + conj(Uold_sub) * conj(Z_N)) * (inv(L_sub))';
    end

    % =========================================================================
    % 5. EVALUACIÓN DE LA ENERGÍA DEL NUEVO ESTADO CANDIDATO
    % =========================================================================
    rho = conj(V) * V.';
    kappa =conj(V) * U.';
    N_final = real(trace(rho));

    % Recalcular Gamma y Delta en el nuevo punto
    Gamma = zeros(m,m);
    for a=1:m
        for c=1:m

            sum_val = 0;

            for d=1:m
                for b=1:m

                    sum_val = sum_val + v_barra(a,b,c,d)*rho(d,b);

                end
            end

            Gamma(a,c) = sum_val;

        end
    end

    % Delta actual
    Delta = zeros(m,m);
    for a=1:m
    for b=1:m
        sum_val = 0;

            for c=1:m
            for d=1:m

                sum_val = sum_val + v_barra(a,b,c,d)*kappa(c,d);

            end
            end

            Delta(a,b) = 0.5*sum_val;

    end
    end
    

    E_actual = trace(t * rho) + 0.5 * trace(Gamma * rho) - 0.5 * trace(Delta * conj(kappa)); 
   % Monitorear el progreso en consola
    fprintf('%8d | %12.6f | %15.4f | %10.2e | %d\n', ...
            iter_E, E_actual, N_final, E_old-E_actual, sub_iter);

    %Verificación de unitariedad Bogoliubov
    tol_unitarity = 1e-7;
    condic1 = norm(U'*U + V'*V - eye(m));
    condic2 = norm(U.'*V + V.'*U);

    if condic1 > tol_unitarity || condic2 > tol_unitarity
        warning('Iteración %d: ¡Pérdida de unitariedad detectada!', k);
        fprintf('  Error Condición 1: %.2e\n', condic1);
        fprintf('  Error Condición 2: %.2e\n', condic2);
    end       

    % Criterio de parada absoluto (convergencia real)
    if abs(E_old-E_actual) < tol_E && abs(N_final - N0) < tol_N
        fprintf('-----------------------------------------------------------------------------------\n');
        fprintf('¡Convergencia HFB absoluta alcanzada en la iteración %d con E = %.6f!\n', iter_E, E_actual);
        break;

    end
    historial_E = [historial_E, real(E_actual)];
    historial_iter = [historial_iter, iter_E];
end