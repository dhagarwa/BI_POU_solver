function [] = helmholtz_sphere()
m =8;
n = 8;
R = 1;
k = 0; %k = 0 laplace
n_patch = 6;
gma = 1;

%Coarse Grid on u-v one patch only
[x1,y1]=ndgrid(pi*(1:(m-2))/(m-1),pi*(1:n-2)/(n-1)); %x1 and y1 both are m-2 x n-2 matrices
nodes =[x1(:),y1(:)]; % N by 2 matrix listing x,y coordinates of all N=mn nodespi
h1 = pi/ (m-1);
h2 = pi/ (n-1);
N = (m-2)*(n-2); %no of nodes on one patch in u-v coordinates


%Singularity resolution grid
hf1 = h1^(3/2);
hf2 = h2^(3/2);
[x2, y2] = ndgrid(-h1: hf1:h1, 0:hf2:pi); %x2 is 
ms = size(x2, 1);
ns = size(x2, 2);
nodes_sing  = [x2(:), y2(:)];
N_sing = size(nodes_sing, 1)

%Fill up load vector
b = [];
for pp = 1:n_patch %pp denotes patch number
    bb = boundary_val(pp, nodes);
    b = [b;bb];
end

%FFT interpolation grid parameters
 uf = 16;  %upsampling factor
 mf = uf*(m-2);
 nf = uf*(n-2);
[xf, yf] = meshgrid(pi*(0:mf-1)/ (mf), pi*(0:nf-1)/ (nf)); %fft upsampled grid in meshgrid format

%Start GMRES
%sigma = gmres(@helmholtz_matvec, b, 3, 1e-5, 20);


%Type the values to file
% sigma_real = real(sigma);
% sigma_imag = imag(sigma);
% adjacent = [sigma_real sigma_imag];
% fileID = fopen('sigma_sphere_801.bin', 'w');
% fwrite(fileID,adjacent,'double');
% fclose(fileID);


% %Integral testing part
% r00 = chart(1, nodes(1, :)); %target node
% % sigma_p = ones(N,1); %sigmas on one patch
% % sigma = [];
% % for pp = 1:n_patch
% %    sigma = [sigma;sigma_p];
% % end
% fileID = fopen('sigma_sphere_800.bin');
% sigma_real = fread(fileID, [N*n_patch, 1], 'double');
% sigma_imag = fread(fileID, [N*n_patch, 1], 'double');
% fclose(fileID);
% sigma = complex(sigma_real, sigma_imag);
% %FFT interpolation preprocessing
% %sigma_pmat = reshape(sigma_p, [m-2, n-2])'; %density in matrix form
%  uf = 16;  %upsampling factor
%  mf = uf*(m-2);
%  nf = uf*(n-2);
% %sigma_pmat = interpft2(sigma_pmat, nf, mf); % density values fft upsampled
%  [xf, yf] = meshgrid(pi*(0:mf-1)/ (mf), pi*(0:nf-1)/ (nf)); %fft upsampled grid in meshgrid format
%  integrate_patches(r00, sigma, nodes, xf, yf, 1)
% 
% 
% %integral testing part ends

%With sigma calculate the solution at some point

fileID = fopen('sigma_sphere_801.bin');
sigma_read_real = fread(fileID, [N*n_patch, 1], 'double');
sigma_read_imag = fread(fileID, [N*n_patch, 1], 'double');
fclose(fileID);
sigma_read = complex(sigma_read_real, sigma_read_imag);



dist = 1.02:0.05:30;
sol = zeros(size(dist, 2), 1);

for kk = 1:size(dist, 2)
    r00 = [0;dist(kk);0]; %target node
    if (dist(kk)-1) < h1
        sol(kk) =  integrate_patches(r00, sigma_read, nodes,  xf, yf, 3);
    elseif (dist(kk)-1) < 5*sqrt(h1)
        sol(kk) =  integrate_patches(r00, sigma_read, nodes,  xf, yf, 2);
    else
        sol(kk) =  integrate_patches(r00, sigma_read, nodes,  xf, yf, 0);
    end
end
%sol(1) = sol(1) + 0.5*sigma_read(N/2);
real_sol = real(sol);
imag_sol = imag(sol);
real_true_sol =  cos(k*dist)./dist;
imag_true_sol = sin(k*dist)./dist;
plot(dist, real_sol, '-o', dist, real_true_sol);

legend('approximate', 'true');
max(abs(real_sol' - real_true_sol))
%(abs(imag_sol - imag_true_sol))
%true_sol./real_sol;


%GMRES matvec function
    function Ax = helmholtz_matvec(x)
        Ax = zeros(N*n_patch, 1);


        for p = 1:n_patch
            for ii = 1:N
                r0 = chart(p, nodes(ii, :));
                Ax((p-1)*N + ii) = 0.5*x((p-1)*N + ii)  + integrate_patches(r0, x, nodes, xf, yf, 1 );            %-0.5x(ii) for interior , + for exterior
            end
        end
    end


%charts for sphere, can be replaced by any surface
function r = chart(patch, node)

    if patch==1 
        u = node(1);
        v = node(2);
        r = [R*sin(u)*cos(v); R*sin(u)*sin(v); R*cos(u)];
    elseif patch==2
        u = node(1);
        v = node(2);
        r = [R*sin(u)*cos(v); -R*sin(u)*sin(v); R*cos(u)];
    elseif patch==3
        u = node(1);
        v = node(2);
        r = [R*sin(u)*cos(v); R*cos(u); R*sin(u)*sin(v)];

    elseif patch==4
        u = node(1);
        v = node(2);
        r = [R*sin(u)*cos(v); R*cos(u); -R*sin(u)*sin(v)];

    elseif patch==5
        u = node(1);
        v = node(2);
        r = [R*sin(u)*sin(v); R*sin(u)*cos(v); R*cos(u)];

    
    elseif patch==6
        u = node(1);
        v = node(2);
        r = [-R*sin(u)*sin(v); R*sin(u)*cos(v); R*cos(u)];
    end 
end

%Exterior Normal to torus
function nhat = normal_sphere(r)
        nhat = r/R;
        
end


%Function to calculate the floating partition of unity
function val_pou = fpou(t_node, s_node, flag)
    d = h1;    %Check if d is correct
    t = norm(t_node - s_node)/ d;
    if t >= 1
        val_pou = 0;
    else 
%         disp('hello');
%         t
        val_pou = exp((2*exp(-1/t))/ (t-1));
    end
    if flag==0
        %disp('hello');
        val_pou =0;
        
    end
end

%Functions to calculate static partition of unity
function val = spou(patch, r)
    patch; r;
    [v, u, rho] = cart2sph(r(1), r(2), r(3));
    %u = pi/2 - u;
    %u; v;
    if patch==1
        r0 = [0;R;0];
    elseif patch ==2    
        r0 = [0;-R;0];
    elseif patch==3
        r0 = [0;0;R];
    elseif patch==4
        r0 = [0;0;-R];
    elseif patch==5
        r0 = [R;0;0];
    elseif patch==6
        r0 = [-R;0;0];
        
    end
    r0;
    [v0, u0, rho0] = cart2sph(r0(1), r0(2), r0(3));
    %u0 = pi/2-u0;
    %u0; v0;
    d = (5/12) * pi*R;
    t = greatCircleDistance(u0, v0, u, v, R)/d;
    if t>=1
        val =0;
    else 
        val = exp((2*exp(-1/t))/ (t-1));
    end
    val;
    
end


%Function to calculate the kernels of doule layer potentials and single
%layer
function val = G1(r0, r1)
    
    val = 1/(4*pi) * cos(k*norm(r0-r1))/ norm(r0-r1);
    %val = 1;
end


function val = G2(r0, r1)

    val = 1i/(4*pi) * sin(k*norm(r0-r1))/ norm(r0-r1);
    %val = 1;
    
end

function val = G3(r0, r1, u, v)
    
    val = 1i/(4*pi) * (sin(k*norm(r0-r1))/ norm(r0-r1)  - k*cos(k*norm(r0-r1))) * (sum((r0-r1).* normal_sphere(r1))/norm(r0-r1)^2);
    %val = 1;

end

function val = G4(r0, r1, u, v)

    val = 1/(4*pi) * (cos(k*norm(r0-r1))/ norm(r0-r1)) * (sum((r0-r1).* normal_sphere(r1))/ norm(r0-r1)^2);
    %val = 1;

end
function val = G5(r0, r1, u, v)

    val = 1/(4*pi) * (k*sin(k*norm(r0-r1))) * (sum((r0-r1).* normal_sphere(r1))/ norm(r0-r1)^2);
    %val = 1;

end

%Function to calculate the green's function
function val = phi(r0, r1)
    g1 = G1(r0, r1);
    g2 = G2(r0, r1);
    val = g1 + g2;
    %val = 100;

end

%Function to calculate the normal derivative of green's function 
function val = dphi(r0, r1, u, v)
    g3 = G3 (r0, r1, u, v);
    g4 = G4(r0, r1, u, v);
    g5 = G5(r0, r1, u, v);
    val = g3 + g4 + g5;
    %val = 100;
    %val = normal_torus(c, a, u, v);
    %val = val(1);
end

%Jacobian for the parameterisation of torus
function val_jac = jacobian(node)
    u = node(1);
    v = node(2);
    val_jac = abs(R^2*sin(u));

end

function val = in_patch(patch, r) %Function to find out if r is in patch no 'patch'
    if patch==1
        r0 = [0;R;0];
    elseif patch ==2    
        r0 = [0;-R;0];
    elseif patch==3
        r0 = [0;0;R];
    elseif patch==4
        r0 = [0;0;-R];
    elseif patch==5
        r0 = [R;0;0];
    elseif patch==6
        r0 = [-R;0;0];
        
    end
        
    [v0, u0, rho0] = cart2sph(r0(1), r0(2), r0(3));
    u0 = pi/2 - u0;
    [v, u, rho] = cart2sph(r(1), r(2), r(3));
    u = pi/2 - u;
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
        r0 = [0;0;R];
        u = angle(r, r0);
        v = angle(r - norm(r)*cos(u)*[0;0;1], [1;0;0]);
    elseif patch ==2
        r0 = [0;0;R];
        u = angle(r0, r);
        v = angle(r - norm(r)*cos(u)*[0;0;1], [1;0;0]);
    elseif patch==3
        r0 = [0;R;0];
        u = angle(r0, r);
        v = angle(r - norm(r)*cos(u)*[0;1;0], [1;0;0]);
    elseif patch==4
        r0 = [0;R;0];
        u = angle(r0, r);
        v = angle(r-norm(r)*cos(u)*[0;1;0], [1;0;0]);
    elseif patch==5
        r0 = [0;0;R];
        u = angle(r0, r);
        v = angle(r - norm(r)*cos(u)*[0;0;1], [0;1;0]);
    elseif patch==6
        r0 = [0;0;R];
        u = angle(r0, r);
        v = angle(r - norm(r)*cos(u)*[0;0;1], [0;1;0]);
        
        
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
        L = 8; % no. of nodes for 1D interpolation
        near_vals = zeros(L+1, 1); % array to store near values
        near_vals(1) = boundary_val(no_patch, no_node);
        sep = zeros(L+1, 1); %separation on x axis between nodes for 1D interpolation
        beta = 1.1;
        for l=2:L+1
            rr = no_r + l*(r0- no_r)*beta*h1/(2*norm(r0-no_r));
            near_val = integrate_patches(rr, x, nodes, xf, yf, 2);
            near_vals(l) = near_val;
            sep(l) = l*beta*h1/2;
        end
        
        sep_q = norm(r0 - no_r);
        val = interp1(sep, near_vals, sep_q, 'spline');

    else
    
        for p=1:n_patch %Iterate over patches
            %p
            if flag ==1  %if r00 is on the surface
                %disp('on surface');
                if in_patch(p, r0) ==1
                    %disp('in patch '), p
                    target_node = get_uv(p, r0);
                    x_p = x((p-1)*N+1: p*N);
                    x_pmat = reshape(x_p, [m-2, n-2])'; %density in matrix form
                    x_pmat = interpft2(x_pmat, nf, mf); % density values fft upsampled
                    val = val + nonsingular_integrate(r0, p, target_node, x_p, nodes, 1);
                    val = val + singular_integrate(p, target_node, x_p, nodes, x_pmat,xf, yf, 1);
                else
                    x_p = x((p-1)*N+1: p*N);
                    val = val + nonsingular_integrate(r0, p, [0, 0], x_p, nodes, 0); %with flag 0 => no fpou , target-node doesn't matter
                    
                end
            elseif flag ==0 %r0 is not on surface , same as well separated integration in evals
                x_p = x((p-1)*N+1: p*N);
                val = val + nonsingular_integrate(r0, p, [0, 0], x_p, nodes, 0); %with flag 0 => no fpou , target-node doesn't matter
                
            elseif flag==2  %Intermediate region integration with fft upsampled values
                x_p = x((p-1)*N+1:p*N);
                x_pmat = reshape(x_p, [m-2, n-2])'; %density in matrix form
                x_pmat = interpft2(x_pmat, nf, mf); % density values fft upsampled
                val = val + inter_integrate(r0, p, x_pmat);
            
            
            
            end
        
        end
    
    end
end






%Non singular integration function
function val = nonsingular_integrate( r0, patch, target_node, x, nodes, flag) %r0, ii(index for u-v corrdinate in nodes array) are both target node coordinates
    %x is density values for one patch only here
    
    val = 0;
    %N = size(nodes, 1);
    %target_node = nodes(ii, :);
    %r0 = chart(target_node);
    for jj= 1:N
        %disp('Iteration over source node '), jj
       source_node = nodes(jj, :);
       if chart(patch, source_node) == r0
           continue   % change here
       end
       u = source_node(1);
       v = source_node(2);
       r1 = chart(patch, source_node);
       wk = spou(patch, r1);
       qk = 0;
       for p= 1:n_patch
           qk = qk+ spou(p, r1);
       end
       
       
%        fpou(target_node, source_node, flag)
%        phi(r0, r1)
%        dphi(r0, r1, u, v)
%        qk
%        jacobian(source_node)
       val = val + h1*h2*(wk/qk)*(-1i*gma*(1 - fpou(target_node, source_node, flag))*phi(r0, r1) * jacobian(source_node) * x(jj) ) + h1*h2*(wk/qk)*(1 - fpou(target_node, source_node, flag))* dphi(r0, r1, u, v) * jacobian(source_node)*x(jj);
        
    end

end


function val = inter_integrate(r0, patch, sigma_pmat)
   val = 0;
   %sigma_uf = sigma_mat(:); % Upsampled sigma by fft
   %N_uf = size(sigma_uf, 1);  % no. of grid points in upsampled grid of sigma_uf 
   h1_uf = pi/mf;
   h2_uf = pi/nf;
   for jj=1:mf
       for ll=1:nf 
      u = xf(1, jj);
      v = yf(ll, 1);
      source_node= [u, v];
      r1 = chart(patch, source_node);
      wk = spou(patch, r1);
      qk = 0;
      for p= 1:n_patch
          qk = qk+ spou(p, r1);
      end
      val = val + h1_uf*h2_uf*(wk/qk)*(-1i*gma)*phi(r0, r1) * jacobian(source_node) * sigma_pmat(ll, jj)  + h1_uf*h2_uf* (wk/qk)*dphi(r0, r1, u, v) * jacobian(source_node)*sigma_pmat(ll,jj); 
       end
   end
    
    
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

% function xval = get_x(source_node, x)
%     u = source_node(1);
%     v = source_node(2);
% %     X = reshape(x, [m-2, n-2])';
% %     uf = 16;
% %     X = interpft2(X, uf*(m-2), uf*(n-2));
% %     mf = uf*(m-2);
% %     nf = uf*(n-2);
% %     [xf, yf] = meshgrid(2*pi*(1:mf)/ (mf+1), 2*pi*(1:nf)/ (nf+1));
%     %size(X ), size(xf), size(yf)
%     
%     xval = interp2(xf, yf, X, u, v, 'spline');
% end


    function index = get_ij(jj)
       
        j = rem(jj, ms);
        i = floor(jj/ ms);
        if j > 0
            i = i+1;
            
        end
        if j==0
           j =  ms; 
        end
        index = [i, j];
        
    end

function val = singular_integrate(patch, target_node, x, nodes, X, xf, yf, flag)
    %%disp('in singular integrate with '), x,ii
    val = 0;
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
    for jj =1 :N_sing
       index = get_ij(jj);  %Using singular node no, get the index of the value we are loking for in interpolated matrix
       source_node_sing = nodes_sing(jj, :);
       rho = source_node_sing(1);
       theta = source_node_sing(2);
       source_node = sing_chart(source_node_sing, target_node);
       u = source_node(1);
       v = source_node(2);
       r1 = chart(patch, source_node);
       wk = spou(patch, r1);
       qk = 0;
       for p= 1:n_patch
           qk = qk+ spou(p, r1);
       end
       
       
       %wk=1;qk=1;
       val = val +  hf1*hf2*(wk/qk)*(-1i*gma*(fpou(target_node, source_node, flag))*phi(r0, r1) * jacobian(source_node) * xval(index(1), index(2)) * abs(rho) ) ;
       val = val + hf1*hf2*(wk/qk)*(fpou(target_node, source_node, flag))* dphi(r0, r1, u, v) * jacobian(source_node)*xval(index(1), index(2))*abs(rho);
        
    end

end



% Function that gives the given boundary Dirchlet data
function val = boundary_val(patch, nodes)
    u = nodes(:, 1);
    v = nodes(:, 2);
    val = ones(size(u, 1), 1) *(exp(1i*k*R)/ R) ;
end


end
    



