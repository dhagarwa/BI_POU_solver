function [] =  helmholtz_evals()

%-----------------------Setting up the variables-------------------------
m =25;
n =25;
c = 2; 
a = 1;
k = 1; %k = 0 laplace
gma = 1;

%------------------------------making the coarse grid----------------------
[x1,y1]=ndgrid(2*pi*(1:(m-2))/(m-1),2*pi*(1:n-2)/(n-1)); %x1 and y1 both are m-2 x n-2 matrices
nodes =[x1(:),y1(:)]; % N by 2 matrix listing x,y coordinates of all N=mn nodespi
h1 = 2*pi/ (m-1);
h2 = 2*pi/ (n-1);
N = (m-2)*(n-2);


%------------------Singularity resolution grid----------------------------
hf1 = h1^(3/2)/2;
hf2 = h2^(3/2)/2;
[x2, y2] = ndgrid(-h1: hf1:h1, 0:hf2:pi); %x2 is 
ms = size(x2, 1);
ns = size(x2, 2);
nodes_sing  = [x2(:), y2(:)];
N_sing = size(nodes_sing, 1);

%--------------------------------Reading the solution sigma--------------- 
fileID = fopen('sigma.txt','r');
formatSpec = '%f';
sigma_read = fscanf(fileID,formatSpec);
fclose(fileID); 
% --------------------------Reading done----------------------------------

% ----------------Code to calculate the solution at points ---------------
dist = 3:0.02:100;
sol = zeros(size(dist, 2), 1);
sigma_mat = reshape(sigma_read, [m-2, n-2])'; %density in matrix form
uf_read = 16;  %upsampling factor
mf_read = uf_read*(m-2);
nf_read = uf_read*(n-2);
sigma_mat = interpft2(sigma_mat, nf_read, mf_read); %nf_read x mf_read % density values fft upsampled 16 times, for inter_integrate h^(3/2) recommended but well i am lazy
[xf_read, yf_read] = meshgrid(2*pi*(0:mf_read-1)/ (mf_read-1), 2*pi*(0:nf_read-1)/ (nf_read-1)); %fft upsampled grid in meshgrid format
for kk = 1:size(dist, 2)
    r00 = [0;dist(kk);0]; %target node
    if (dist(kk) - 3) <0.1
        sol(kk) = near_integrate(r00, sigma_mat, nodes);
    elseif (dist(kk) - 3) < 2
        sol(kk) =inter_integrate(r00,sigma_mat, nodes);
    else
        sol(kk) = nonsingular_integrate(r00, 1, sigma_read, nodes, 0);
    end
end
real_sol = real(sol);
plot(dist-3, real_sol, '-o');

%--------------------------Calculation and plotting ends------------------

%--------------------------Essential functions----------------------------
%chart for torus, can be replaced by any surface
function r = chart(node)

    u = node(1);
    v = node(2);
    r = [(c + a*cos(v))*cos(u); (c + a*cos(v))*sin(u); a*sin(v)]; 

end

%Exterior Normal to torus
    function nhat = normal_torus(c, a, u, v)
        
        rout = c;
        rin = a;
        iangle = v;
        jangle = u;
%         vx = cos(jangle)*(rout+cos(iangle)*rin);
%         vy = sin(jangle)*(rout+cos(iangle)*rin);
%         vz = sin(iangle)*rin;
        
        tx = -sin(jangle);
        ty = cos(jangle);
        tz = 0;
        
        sx = cos(jangle)*(-sin(iangle));
        sy = sin(jangle)*(-sin(iangle));
        sz = cos(iangle);
        
        nx = ty*sz - tz*sy;
        ny = tz*sx - tx*sz;
        nz = tx*sy - ty*sx;
        
        length = sqrt(nx*nx + ny*ny + nz*nz);
        nx = -nx/ length;
        ny = -ny / length;
        nz = -nz/ length;
        nhat = [nx; ny; nz];
        
    end


%Function to calculate the floating partition of unity
function val_pou = fpou(t_node, s_node, flag)
    d = h1;    %Check if d is correct
    t = norm(t_node - s_node)/ d;
    if t >= 1
        val_pou = 0;
    else 
        val_pou = exp((2*exp(-1/t))/ (t-1));
    end
    if flag==0
        val_pou =0;
        
    end
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
    
    val = 1i/(4*pi) * (sin(k*norm(r0-r1))/ norm(r0-r1)^3  - k*cos(k*norm(r0-r1))/ norm(r0-r1)^2) * sum((r0-r1).* normal_torus(c, a, u, v));
    %val = 1;

end

function val = G4(r0, r1, u, v)

    val = 1/(4*pi) * (cos(k*norm(r0-r1))/ norm(r0-r1)^3) * sum((r0-r1).* normal_torus(c, a, u, v));
    %val = 1;

end
function val = G5(r0, r1, u, v)

    val = 1/(4*pi) * (k*sin(k*norm(r0-r1))/ norm(r0-r1)^2) * sum((r0-r1).* normal_torus(c,a, u, v));
    %val = 1;

end

%Function to calculate the green's function
function val = phi(r0, r1)
    g1 = G1(r0, r1);
    g2 = G2(r0, r1);
    val = g1 + g2;
    %val = 1;

end

%Function to calculate the normal derivative of green's function 
function val = dphi(r0, r1, u, v)
    g3 = G3 (r0, r1, u, v);
    g4 = G4(r0, r1, u, v);
    g5 = G5(r0, r1, u, v);
    val = g3 + g4 + g5;
    %val = normal_torus(c, a, u, v);
    %val = val(1);
end

%Jacobian for the parameterisation of torus
function val_jac = jacobian(node)
    u = node(1);
    v = node(2);
    val_jac = abs(a*(c + a*cos(v)));
end


%-----------------------------------Integration functions start-------------------------------------------------------------- 

% the function for well separated integration(separation > sqrt(h1)): already there just copy it
function val = nonsingular_integrate( r0, ii, sigma, nodes, flag)
    %disp('In non singular integrate with '), x 
    val = 0;
    %N = size(nodes, 1);
    target_node = nodes(ii, :);
    %r0 = chart(target_node);
    for jj= 1:N
        %disp('Iteration over source node '), jj
       source_node = nodes(jj, :);
       if chart(source_node) == r0
           continue   % change here
       end
       u = source_node(1);
       v = source_node(2);
       r1 = chart(source_node);
       val = val + h1*h2*(-1i*gma*(1 - fpou(target_node, source_node, flag))*phi(r0, r1) * jacobian(source_node) * sigma(jj) ) + h1*h2*(1 - fpou(target_node, source_node, flag))* dphi(r0, r1, u, v) * jacobian(source_node)*sigma(jj);
        
    end

end

%Well separated integration ends
%-------------------------------------------------------------------------

%Intermediate separation integration( h1<dist < sqrt(h1) )
function val = inter_integrate(r0, sigma_mat, nodes)
    
   val = 0;
   sigma_uf = sigma_mat(:); % Upsampled sigma by fft
   N_uf = size(sigma_uf, 1);  % no. of grid points in upsampled grid of sigma_uf 
   h1_uf = 2*pi/mf_read;
   h2_uf = 2*pi/nf_read;
   for jj=1:mf_read
       for ll=1:nf_read 
      u = xf_read(1, jj);
      v = yf_read(ll, 1);
      source_node= [u, v];
      r1 = chart(source_node);
      val = val + h1_uf*h2_uf*(-1i*gma)*phi(r0, r1) * jacobian(source_node) * sigma_mat(ll, jj)  + h1_uf*h2_uf* dphi(r0, r1, u, v) * jacobian(source_node)*sigma_mat(ll,jj); 
       end
   end
    
    
end


%Intermediate integration ends

%----------------------------------------------------------------------------
% nearest integration (dist < h1)
function val = near_integrate(r0, sigma_mat, nodes)
    no_node = get_no(r0); %get near orthogonal node parametrised
    no_r = chart(no_node);
    L = 8; % no. of nodes for 1D interpolation
    near_vals = zeros(L+1, 1);
    near_vals(1) = boundary_val(no_node);
    sep = zeros(L+1, 1);
    beta = 1.1;
    for l=2:L+1
        rr = no_r + l*(r0- no_r)*beta*h1/(3*norm(r0-no_r));
        near_val = inter_integrate(rr, sigma_mat, nodes)
        near_vals(l) = near_val;
        sep(l) = l*beta*h1/3;
    end
    
    sep_q = norm(r0 - no_r);
    val = interp1(real(sep), real(near_vals), sep_q, 'spline');
end

%nearest integration ends



%singular integration


function val = singular_integrate(ii, x, nodes, X, xf, yf, flag)
    %%disp('in singular integrate with '), x,ii
    val = 0;
%     hf1 = h1^(3/2)/3;
%     hf2 = h2^(3/2)/3;
%     [x2, y2] = ndgrid(-h1+hf1: hf1:h1-hf1, 0+hf2:hf2:pi-hf2);
%     nodes_sing  = [x2(:), y2(:)];
%     N_sing = size(nodes_sing, 1);
    
    target_node = nodes(ii, :);
    r0 = chart(target_node);
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
       r1 = chart(source_node);
       %xval = get_x(source_node, x);
       val = val +  hf1*hf2*(-1i*gma*(fpou(target_node, source_node, flag))*phi(r0, r1) * jacobian(source_node) * xval(index(1), index(2)) * abs(rho) ) ;
       val = val + hf1*hf2*(fpou(target_node, source_node, flag))* dphi(r0, r1, u, v) * jacobian(source_node)*xval(index(1), index(2))*abs(rho);
        
    end

end
%----------------------Some helper functions------------------------------

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


    
    
 %get near orthogonal node   
function no_node = get_no(r)
    no_scale = @(xx) dot(r - chart(xx), normal_torus(c, a, xx(1), xx(2)))/ norm(r - chart(xx));
    ini_val = [pi/2, 0.1];
    [no_node, fval] = fminsearch(no_scale , ini_val);
    
        
        
end


function val = boundary_val(nodes)
    u = nodes(:, 1);
    v = nodes(:, 2);
    val = ones(size(u, 1), 1) ;
end


end


