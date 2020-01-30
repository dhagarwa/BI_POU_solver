function [] = helmholtz_sphereVec(m, n, k, gma)
%m =30;
%n = 30;
R = 1;
%k = 1; %k = 0 laplace
n_patch = 6;
%gma = 1;

%Coarse Grid on u-v one patch only
[x1,y1]=ndgrid(pi*(1:(m-2))/(m-1),pi*(1:n-2)/(n-1)); %x1 and y1 both are m-2 x n-2 matrices
nodes =[x1(:),y1(:)]; % N by 2 matrix listing x,y coordinates of all N=mn nodespi
h1 = pi/ (m-1);
h2 = pi/ (n-1);
N = (m-2)*(n-2); %no of nodes on one patch in u-v coordinates


%Singularity resolution grid
hf1 = h1^(3/2);
hf2 = h2^(1/1);
rho_sing = 3*sqrt(h1);
%hf1 =  rho_sing/128;
%hf2 = pi/128;
[x2, y2] = ndgrid(-rho_sing :hf1:rho_sing-hf1, 0:hf2:2*pi-hf2); %x2 is 
ms = size(x2, 1)
ns = size(x2, 2)
nodes_sing  = [x2(:), y2(:)];
N_sing = size(nodes_sing, 1)

%Fill up load vector
b = [];
for pp = 1:n_patch %pp denotes patch number
    bb = boundary_val(pp, nodes);
    b = [b;bb];
end

%FFT interpolation grid parameters
 uf = 4;  %upsampling factor
 mf = uf*(m-2);
 nf = uf*(n-2);
[xf, yf] = meshgrid(pi*(0:mf-1)/ (mf), pi*(0:nf-1)/ (nf)); %fft upsampled grid in meshgrid format

%Start GMRES
%sigma = gmres(@helmholtz_matvec, b, 3, 1e-10, 20);

s1 = 'sigma_vec_int_';
s2 = num2str(m);
s3 = num2str(k);
s4 = '.bin';
lbl = strcat(s1, s2, s3, s4);

%Type the values to file
% sigma_real = real(sigma);
% sigma_imag = imag(sigma);
% adjacent = [sigma_real sigma_imag];
% fileID = fopen(lbl, 'w');
% fwrite(fileID,adjacent,'double');
% fclose(fileID);


% %------Integral testing part---------------------------

% sigma_p = ones(N,1); %sigmas on one patch
% sigma = [];
% for pp = 1:n_patch
%    sigma = [sigma;sigma_p];
% end
% 
sigma_read = readSigma(lbl, N, n_patch);
r00 = [0.0, 0.99 , 0]; %target node
uf = 4;  %upsampling factor
mf = uf*(m-2);
nf = uf*(n-2);
[xf, yf] = meshgrid(pi*(0:mf-1)/ (mf), pi*(0:nf-1)/ (nf)); %fft upsampled grid in meshgrid format

area = 1 + r00(2);

err1 = abs(integrate_patches(r00, sigma_read, nodes, xf, yf, 3) - area)

% 
% % 
% % 
% %------------integral testing part ends-------------


% % -------------------Equivalent sources code------------------
% %Calculating values at inward check surface 
% sigma_read = readSigma(lbl, N, n_patch);
% check_points = [];
% equi_points = [];
% R_equi = 0.90;
% R_check = 0.75;
% azi_h = 2*pi/(3*m); polar_h = pi/(3*n);
% azi = 0:azi_h:2*pi-0.1; polar = 0:polar_h:pi;
% for count1=1:length(azi)
%     for count2 = 1:length(polar)
%         check_x = R_check*cos(azi(count1))*sin(polar(count2));
%         check_y = R_check*sin(azi(count1))* sin(polar(count2));
%         check_z = R_check*cos(polar(count2));
%         check_points = [check_points; check_x, check_y, check_z ];
%         equi_x = R_equi*cos(azi(count1))*sin(polar(count2));
%         equi_y = R_equi*sin(azi(count1))* sin(polar(count2));
%         equi_z = R_equi*cos(polar(count2));
%         equi_points = [equi_points; equi_x, equi_y, equi_z ];
%     end
% end
% 
% 
% values_check_surface = zeros(size(check_points, 1), 1);
% for count = 1:size(check_points, 1)
%     values_check_surface(count) = integrate_patches(check_points(count, :), sigma_read, nodes, xf, yf, 2);
% end
% 
% %values_check_surface
% 
% %Use GMRES to solve for "w" : equivalent weights at check points 
% weights = gmres(@equi_source_matvec, values_check_surface, 5, 1e-10, 30);
% trg_point = [0.2, 0.2, 0];
% disp('Error:');
% abs(sum(dphi(trg_point, equi_points).*weights) - (1+trg_point(2)))
 % % ---------------------Equivalent sources code ends-------------------

% % -------------Calculating solution, error and graph -----------------
% %With sigma calculate the solution at some point
% 
% fileID = fopen(lbl);
% sigma_read_real = fread(fileID, [N*n_patch, 1], 'double');
% sigma_read_imag = fread(fileID, [N*n_patch, 1], 'double');
% fclose(fileID);
% sigma_read = complex(sigma_read_real, sigma_read_imag);
% 
% 
% [x_err, y_err, z_err] = meshgrid(0:0, 0:0.05:0.5, 0:0 ); %error grid x and y
% weight = (1/10)^3;
% nodes_err = [x_err(:), y_err(:), z_err(:)];
% N_err = size(nodes_err, 1) ;
% total_err = 0;
% total_err_inf = 0;
% solution = 0;
% for iii = 1:N_err
%    node_err = nodes_err(iii, :);
%    dist = norm(node_err);
%    r00 = [0, dist, 0];
%    if abs(dist-R) < h1
%        sol = integrate_patches(r00, sigma_read, nodes, xf, yf, 3);
%    elseif abs(dist-R) < sqrt(h1)
%        sol = integrate_patches(r00, sigma_read, nodes, xf, yf, 0);
%    else
%        sol = integrate_patches(r00, sigma_read, nodes, xf, yf, 0);
%    end
%    
%    %true_sol = cos(k*dist)/dist + 1i*sin(k*dist)/dist;
%    true_sol =1+dist;
%    total_err = total_err + (norm(sol-true_sol))^2 * weight;
%    total_err_inf = max(total_err_inf, abs(sol-true_sol));
%    solution = solution + norm(true_sol)^2 *weight;
%    
% end
% m, k, gma, h1
% %total_err = total_err^(0.5) / solution^(0.5)
% total_err_inf
% dist = 0:0.05:1;
% sol = zeros(size(dist, 2), 1);
% 
% for kk = 1:size(dist, 2)
%     r00 = [0,dist(kk),0]; %target node
%     if abs(dist(kk)-R) < h1
%         sol(kk) =  integrate_patches(r00, sigma_read, nodes,  xf, yf, 3);
%     elseif abs(dist(kk)-R) < sqrt(h1)
%         sol(kk) =  integrate_patches(r00, sigma_read, nodes,  xf, yf, 0);
%     else
%         sol(kk) =  integrate_patches(r00, sigma_read, nodes,  xf, yf, 0);
%     end
% end
% 
% real_sol = real(sol);
% imag_sol = imag(sol);
% real_true_sol = 1+dist;
% %real_true_sol =  cos(k*dist)./dist;
% %imag_true_sol = sin(k*dist)./dist;
% plot(dist, real_sol, '-o', dist, real_true_sol);
% 
% legend('approximate', 'true');
% xlabel('r');
% ylabel('Realpart of solution');
% % -------------Calculating solution and graph ends ---------

%GMRES matvec function
    function Ax = helmholtz_matvec(x)
        Ax = zeros(N*n_patch, 1);


        for p = 1:n_patch
            for ii = 1:N
                r0 = chart(p, nodes(ii, :));
                Ax((p-1)*N + ii) = -0.5*x((p-1)*N + ii)  + integrate_patches(r0, x, nodes, xf, yf, 1 );            %-0.5x(ii) for interior , + for exterior
            end
        end
    end


% GMRES equisource matevc function
    function Lx = equi_source_matvec(x)
        N_check = size(values_check_surface, 1);
       Lx = zeros(N_check, 1);
       for ii=1:N_check
          Lx(ii) = sum(dphi(check_points(ii, :), equi_points).*x); 
       end
        
    end


%charts for sphere, can be replaced by any surface
function r = chart(patch, node)

    if patch==1 
        u = node(:, 1);
        v = node(:, 2);
        r = [R*sin(u).*cos(v), R*sin(u).*sin(v), R*cos(u)];
    elseif patch==2
        u = node(:,1);
        v = node(:,2);
        r = [R*sin(u).*cos(v), -R*sin(u).*sin(v), R*cos(u)];
    elseif patch==3
        u = node(:,1);
        v = node(:,2);
        r = [R*sin(u).*cos(v), R*cos(u), R*sin(u).*sin(v)];

    elseif patch==4
        u = node(:,1);
        v = node(:,2);
        r = [R*sin(u).*cos(v), R*cos(u), -R*sin(u).*sin(v)];

    elseif patch==5
        u = node(:,1);
        v = node(:,2);
        r = [R*sin(u).*sin(v), R*sin(u).*cos(v), R*cos(u)];

    
    elseif patch==6
        u = node(:,1);
        v = node(:,2);
        r = [-R*sin(u).*sin(v), R*sin(u).*cos(v), R*cos(u)];
    end 
end

%Exterior Normal to torus
function nhat = normal_sphere(r)
        nhat = r/R;
        
end


%Function to calculate the floating partition of unity
function val_pou = fpou(t_node, s_node, flag)
    d = rho_sing;    %Check if d is correct
    t_node = repmat(t_node, size(s_node, 1), 1);
    t = sqrt(sum((t_node - s_node).^2, 2))/ d;
    t(t>=1) = 100;
    t(t ==0) = 200;
    t(t < 1) = exp((2*exp(-1./t(t < 1)))./ (t(t<1)-1));
    t(t==100) = 0;
    t(t==200) = 1;
%     if t >= 1
%         val_pou = 0;
%     else 
% %         disp('hello');
% %         t
%         val_pou = exp((2*exp(-1/t))/ (t-1));
%     end
    if flag==0
        %disp('hello');
        t(:) =0;
        
    end
    %norm(t, inf)
    %t(:) = 1;
    val_pou = t;
end

%Functions to calculate static partition of unity
function val = spou(patch, r)
    [v, u, rho] = cart2sph(r(:,1), r(:,2), r(:,3));
    %u = pi/2 - u;
    %u; v;
    if patch==1
        r0 = [0,R,0];
    elseif patch ==2    
        r0 = [0,-R,0];
    elseif patch==3
        r0 = [0,0,R];
    elseif patch==4
        r0 = [0,0,-R];
    elseif patch==5
        r0 = [R,0,0];
    elseif patch==6
        r0 = [-R,0,0];
        
    end
    [v0, u0, rho0] = cart2sph(r0(1), r0(2), r0(3));
    %u0 = pi/2-u0;
    %u0; v0;
    d = (5/12) * pi*R;
    nn = size(r, 1);
    val = zeros(nn, 1);
    for ii=1:nn
        t = greatCircleDistance(u0, v0, u(ii), v(ii), R)/d;
        if t>=1
            val(ii) =0;
        elseif t ==0
            val(ii) = 1;
        else
            val(ii) = exp((2*exp(-1/t))/ (t-1));
        end
    end
    
    
end




%Function to calculate the green's function
function val = phi(r0, r1)
    r0 = repmat(r0, size(r1, 1), 1);
    nrm = sqrt(sum((r0-r1).^2, 2));
    nrm(nrm==0) = 1;
    g1 = 1/(4*pi) * cos(k*nrm)./ nrm;
    g2 = 1i/(4*pi) * sin(k*nrm)./ nrm;
    val = g1 + g2;
    %val = ones(size(r1, 1), 1);

end

%Function to calculate the normal derivative of green's function 
function val = dphi(r0, r1)
    r0 = repmat(r0, size(r1, 1), 1);
    nrm = sqrt(sum((r0-r1).^2, 2));
    nrm(nrm==0) = 1; 
    g3 = 1i/(4*pi) * (sin(k*nrm)./ nrm  - k*cos(k*nrm)) .* (sum((r0-r1).* (r1/R), 2)./(nrm.^2));
    g4 = 1/(4*pi) * (cos(k*nrm)./ nrm) .* (sum((r0-r1).* (r1/R), 2)./ (nrm.^2));
    g5 = 1/(4*pi) * (k*sin(k*nrm)) .* (sum((r0-r1).* (r1/R), 2)./ (nrm.^2));
    val = g3 + g4 + g5;
    %val = ones(size(r1, 1), 1);
    
end

%Jacobian for the parameterisation of sphere
function val_jac = jacobian(node)
    u = node(:, 1);
    v = node(:,2);
    val_jac = abs(R^2*sin(u));

end

function val = in_patch(patch, r) %Function to find out if r is in patch no 'patch'
    if patch==1
        r0 = [0,R,0];
    elseif patch ==2    
        r0 = [0,-R,0];
    elseif patch==3
        r0 = [0,0,R];
    elseif patch==4
        r0 = [0,0,-R];
    elseif patch==5
        r0 = [R,0,0];
    elseif patch==6
        r0 = [-R,0,0];
        
    end
        
    [v0, u0, rho0] = cart2sph(r0(1), r0(2), r0(3));
    %u0 = pi/2 - u0;
    [v, u, rho] = cart2sph(r(1), r(2), r(3));
    %u = pi/2 - u;
    d = greatCircleDistance(u, v, u0, v0, R);
    if d >= pi*R/2
        val = 0;
    else 
        val = 1;
    end
        
end

%Function to find the angle between two vectors a, b
function val = angle(a, b)
    val = atan2(norm(cross(a,b)), dot(a,b));
        
end

%Function to calculate the parametric coodinates of a cartesian point in a patch given by patch no patch 
function node = get_uv(patch, r)
    if patch==1
        r0 = [0,0,R];
        u = angle(r, r0);
        v = angle(r - norm(r)*cos(u)*[0,0,1], [1,0,0]);
    elseif patch ==2
        r0 = [0,0,R];
        u = angle(r0, r);
        v = angle(r - norm(r)*cos(u)*[0,0,1], [1,0,0]);
    elseif patch==3
        r0 = [0,R,0];
        u = angle(r0, r);
        v = angle(r - norm(r)*cos(u)*[0,1,0], [1,0,0]);
    elseif patch==4
        r0 = [0,R,0];
        u = angle(r0, r);
        v = angle(r-norm(r)*cos(u)*[0,1,0], [1,0,0]);
    elseif patch==5
        r0 = [0,0,R];
        u = angle(r0, r);
        v = angle(r - norm(r)*cos(u)*[0,0,1], [0,1,0]);
    elseif patch==6
        r0 = [0,0,R];
        u = angle(r0, r);
        v = angle(r - norm(r)*cos(u)*[0,0,1], [0,1,0]);
        
        
    end
        
    node = [u, v];
end


%Non singular integrate iterate over patches
function val = integrate_patches(r0, x, nodes, xf, yf, flag) %x is density in all patches, r0 target nod in cartesian
    val =0;
    %flag=0, not on surface, flag=1, on surface nsi+si, flag=2 intermdiate
    %integration, flag=3, nearest integration
    if flag==3 % nearest integration, no need to go into patches
        no_node = [pi/2, pi/2]; %near orthogonal node
        no_patch = 1; %near orthogonal patch no
        no_r = chart(no_patch, no_node);
        L = 7; % no. of nodes for 1D interpolation
        near_vals = zeros(L+1, 1); % array to store near values
        near_vals(1) = boundary_val(no_patch, no_node);
        sep = zeros(L+1, 1); %separation on x axis between nodes for 1D interpolation
        beta = 1.1;
        for l=2:L+1
            rr = no_r + l*(r0- no_r)*beta*h1/(2*norm(r0-no_r));
            near_val = integrate_patches(rr, x, nodes, xf, yf, 0);
            near_vals(l) = near_val;
            sep(l) = l*beta*h1/2;
        end
        
        sep_q = norm(r0 - no_r);
        val = interp1(sep, near_vals, sep_q, 'cubic');

    else
    
        for p=1:6 %Iterate over patches
            %p
            if flag ==1  %if r00 is on the surface
                %disp('on surface');
                if in_patch(p, r0) ==1
                    %disp('in patch '), p
                    target_node = get_uv(p, r0);
                    r01 = chart(p, target_node);
                    if norm(r01 - r0) > 1e-6
                        disp('get uv incorrect');
                    end
                    x_p = x((p-1)*N+1: p*N);
                    r1 = chart(p, nodes);
                    qk = zeros(N, 1);
                    for ppp= 1:n_patch
                        qk = qk+ spou(ppp, r1);
                    end
                    wk  = spou(p, r1);
                    x_p = x_p.*(wk./qk);        
                    x_pmat = reshape(x_p, [m-2, n-2])'; %density in matrix form
                    x_pmat = [zeros(1, m-2);x_pmat];  
                    x_pmat = [zeros(n-1, 1), x_pmat];
                    x_pmat = interpft2(x_pmat, nf, mf); % density values fft upsampled
                    val = val + nonsingular_integrate(r0, p, target_node, x_p, nodes, 1);
                    val = val + singular_integrate(p, target_node, x_p, nodes, x_pmat,xf, yf, 1);
                else
                    x_p = x((p-1)*N+1: p*N);
                    r1 = chart(p, nodes);
                    qk = zeros(N, 1);
                    for ppp= 1:n_patch
                        qk = qk+ spou(ppp, r1);
                    end
                    wk  = spou(p, r1);
                    x_p = x_p.*(wk./qk);       
                    val = val + nonsingular_integrate(r0, p, [0, 0], x_p, nodes, 0); %with flag 0 => no fpou , target-node doesn't matter
                    
                end
            elseif flag ==0 %r0 is not on surface , same as well separated integration in evals
                x_p = x((p-1)*N+1: p*N);                
                r1 = chart(p, nodes);
                qk = zeros(N, 1);
                for ppp= 1:n_patch
                    qk = qk+ spou(ppp, r1);
                end
                wk  = spou(p, r1);
                x_p = x_p.*(wk./qk); 
                val = val + nonsingular_integrate(r0, p, [0, 0], x_p, nodes, 0); %with flag 0 => no fpou , target-node doesn't matter
                
            elseif flag==2  %Intermediate region integration with fft upsampled values
                x_p = x((p-1)*N+1:p*N);
                r1 = chart(p, nodes);
                qk = zeros(N, 1);
                for ppp= 1:n_patch
                    qk = qk+ spou(ppp, r1);
                end
                wk  = spou(p, r1);
                x_p = x_p.*(wk./qk);
                x_pmat = reshape(x_p, [m-2, n-2])'; %density in matrix form
                x_pmat = [zeros(1, m-2);x_pmat];
                x_pmat = [zeros(n-1, 1), x_pmat];
                x_pmat = interpft2(x_pmat, nf, mf); % density values fft upsampled
                val = val + inter_integrate(r0, p, x_pmat);
                
            
            
            end
        
        end
    
    end
end






%Non singular integration function
function val = nonsingular_integrate( r0, patch, target_node, x, nodes, flag) %r0, ii(index for u-v corrdinate in nodes array) are both target node coordinates
    %x is density values for one patch only here
    r1 = chart(patch, nodes);
    %l = (1 - fpou(target_node, nodes, flag));
    Ssigma = (1 - fpou(target_node, nodes, flag)).*phi(r0, r1).*jacobian(nodes).*x;
    Dsigma = (1 - fpou(target_node, nodes, flag)).*dphi(r0, r1).*jacobian(nodes).*x;
    %Dsigma = (1 - fpou(target_node, nodes, flag)).*dphi(r0, r1).*x;
    val = h1*h2 * (-1i *gma)* sum(Ssigma)   + h1*h2*sum(Dsigma);
end


function val = inter_integrate(r0, patch, sigma_pmat)
   
   %sigma_uf = sigma_mat(:); % Upsampled sigma by fft
   %N_uf = size(sigma_uf, 1);  % no. of grid points in upsampled grid of sigma_uf 
   h1_uf = pi/mf;
   h2_uf = pi/nf;
   source_nodes = [xf(:), yf(:)];
   r1 = chart(patch, source_nodes);
   sigma_pmat = sigma_pmat(:);
   Ssigma = phi(r0, r1).*jacobian(source_nodes).*sigma_pmat;
   
   Dsigma = dphi(r0, r1).*jacobian(source_nodes).*sigma_pmat;
   %Dsigma = dphi(r0, r1).*sigma_pmat;
   val = h1_uf*h2_uf*(-1i*gma)*sum(Ssigma) + h1_uf*h2_uf*sum(Dsigma);
   
     
end 
        

function node = sing_chart(source_node, target_node)
    rho = source_node(:,1);
    theta = source_node(:,2);
    u = target_node(1) + rho.*cos(theta);
    v = target_node(2) + rho.*sin(theta);
    node = [u, v];

end

%2D FFT Interpolation on refined grid
function Y=interpft2(X,p,q)
X = interpft(X,p);
Y = interpft(X',q);
Y = Y';

end

function val = singular_integrate(patch, target_node, x, nodes, X, xf, yf, flag)
    %%disp('in singular integrate with '), x,ii
%     hf1 = h1^(3/2)/3;
%     hf2 = h2^(3/2)/3;
%     [x2, y2] = ndgrid(-h1+hf1: hf1:h1-hf1, 0+hf2:hf2:pi-hf2);
%     nodes_sing  = [x2(:), y2(:)];
%     N_sing = size(nodes_sing, 1);
    
    %target_node = nodes(ii, :);
    r0 = chart(patch, target_node);
    source_nodes = sing_chart(nodes_sing,target_node); %all singular resolution nodes in u-v coordinates
    all_u = reshape(source_nodes(:, 1), [ms, ns])';
    all_v = reshape(source_nodes(:, 2), [ms, ns])';
    xval = interp2(xf, yf, X, all_u, all_v, 'spline'); %interpolated density values
    xval = xval';
    xval = xval(:);
    rho = nodes_sing(:, 1);
    r1 = chart(patch, source_nodes);
    Ssigma = fpou(target_node, source_nodes, flag).*phi(r0, r1).*jacobian(source_nodes).*xval.*abs(rho);
    Dsigma = fpou(target_node, source_nodes, flag).*dphi(r0, r1).*jacobian(source_nodes).*xval.*abs(rho);
    %Dsigma = fpou(target_node, source_nodes, flag).*dphi(r0, r1).*xval.*abs(rho);
    val =0.5*( hf1*hf2*(-1i*gma)*sum(Ssigma) + hf1*hf2*sum(Dsigma));
    %val =0.5*hf1*hf2*sum(Dsigma);
    %val_fake =0.5*( hf1*hf2*(-1i*gma)*sum(Ssigma) + hf1*hf2*sum(Dsigma_fake));
    %err_sing = abs(val_fake - pi*(h1^2))

end



% Function that gives the given boundary Dirchlet data
function val = boundary_val(patch, nodes)
    u = nodes(:, 1);
    v = nodes(:, 2);
    %val = ones(size(u, 1), 1) *(exp(1i*k*R)/ R) ;
    val = zeros(size(nodes, 1), 1);
    for ii=1:size(nodes, 1)
        r = chart(patch, nodes(ii, :));
        val (ii) = 1+r(2);
    end
end


end
    



