function [] = helmholtz2()

m =25;
n =25;
c = 2; 
a = 1;
k = 1; %k = 0 laplace
gma = 1;


[x1,y1]=ndgrid(2*pi*(1:(m-2))/(m-1),2*pi*(1:n-2)/(n-1)); %x1 and y1 both are m-2 x n-2 matrices
nodes =[x1(:),y1(:)]; % N by 2 matrix listing x,y coordinates of all N=mn nodespi
h1 = 2*pi/ (m-1);
h2 = 2*pi/ (n-1);
N = (m-2)*(n-2);


%Singularity resolution grid
hf1 = h1^(3/2)/2;
hf2 = h2^(3/2)/2;
[x2, y2] = ndgrid(-h1: hf1:h1, 0:hf2:pi); %x2 is 
ms = size(x2, 1);
ns = size(x2, 2);
nodes_sing  = [x2(:), y2(:)];
N_sing = size(nodes_sing, 1)

%Fill up load vector 
b = boundary_val(nodes);


%Start GMRES
%sigma = gmres(@helmholtz_matvec, b, 3, 1e-5, 20);


%Type the values to file
% fileID = fopen('sigma.txt','w');
% fprintf(fileID,'%4.4f\n',sigma);
% fclose(fileID);

% %Integral testing part
% r00 = [10; 0; 0]; %target node
% ii = 1;
% x = ones(N,1);
% 
% %FFT interpolation preprocessing
% X = reshape(x, [m-2, n-2])'; %density in matrix form
% uf = 16;  %upsampling factor
% mf = uf*(m-2);
% nf = uf*(n-2);
% X = interpft2(X, nf, mf); % density values fft upsampled
% [xf, yf] = meshgrid(2*pi*(0:mf-1)/ (mf-1), 2*pi*(0:nf-1)/ (nf-1)); %fft upsampled grid in meshgrid format
% nsi = nonsingular_integrate(r00, x, nodes, 1)
% si = singular_integrate(ii, x, nodes, X, xf, yf,1)
% nsi + si
% %integral testing part ends

%With sigma calculate the solution at some point

fileID = fopen('sigma.txt','r');
formatSpec = '%f';
sigma_read = fscanf(fileID,formatSpec);
fclose(fileID); 

dist = -3.0:-0.1:-100;
sol = zeros(size(dist, 2), 1);
        X_read = reshape(sigma_read, [m-2, n-2])'; %density in matrix form
        uf_read = 16;  %upsampling factor
        mf_read = uf_read*(m-2);
        nf_read = uf_read*(n-2);
        X_read = interpft2(X_read, nf_read, mf_read); % density values fft upsampled
        [xf_read, yf_read] = meshgrid(2*pi*(0:mf_read-1)/ (mf_read-1), 2*pi*(0:nf_read-1)/ (nf_read-1)); %fft upsampled grid in meshgrid format
for kk = 1:size(dist, 2)
    r00 = [0;dist(kk);0]; %target node
    if kk==1
        sol(kk) = sol(kk) + singular_integrate(1, sigma_read, nodes, X_read, xf_read, yf_read, 1);
        sol(kk) = sol(kk)+ nonsingular_integrate(r00, 1, sigma_read, nodes, 1);
    else
        sol(kk) =nonsingular_integrate(r00, 1, sigma_read, nodes, 0);
    end
end
real_sol = real(sol);
plot(dist, real_sol);


%GMRES matvec function
    function Ax = helmholtz_matvec(x)
        Ax = zeros(N, 1);
        %FFT interpolation preprocessing
        X = reshape(x, [m-2, n-2])'; %density in matrix form
        uf = 16;  %upsampling factor
        mf = uf*(m-2);
        nf = uf*(n-2);
        X = interpft2(X, nf, mf); % density values fft upsampled
        [xf, yf] = meshgrid(2*pi*(0:mf-1)/ (mf-1), 2*pi*(0:nf-1)/ (nf-1)); %fft upsampled grid in meshgrid format
        for ii = 1:N
            r0 = chart(nodes(ii, :));
            si = singular_integrate(ii, x, nodes, X, xf, yf, 1);
            nsi = nonsingular_integrate(r0, ii, x, nodes, 1);
           Ax(ii) = 0.5*x(ii)  + nsi +si;            %-0.5x(ii) for interior , + for exterior
        end
    end

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
%         disp('hello');
%         t
        val_pou = exp((2*exp(-1/t))/ (t-1));
    end
    if flag==0
        %disp('hello');
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

%Non singular integration function
function val = nonsingular_integrate( r0, ii, x, nodes, flag)
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
       val = val + h1*h2*(-1i*gma*(1 - fpou(target_node, source_node, flag))*phi(r0, r1) * jacobian(source_node) * x(jj) ) + h1*h2*(1 - fpou(target_node, source_node, flag))* dphi(r0, r1, u, v) * jacobian(source_node)*x(jj);
        
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



% Function that gives the given boundary Dirchlet data
function val = boundary_val(nodes)
    u = nodes(:, 1);
    v = nodes(:, 2);
    val = ones(size(u, 1), 1) ;
end


end