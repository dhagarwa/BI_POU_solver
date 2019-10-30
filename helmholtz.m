function [] = helmholtz()

% Solving an exterior Dirichlet problem for Helmholtz equation
% We use the integral equations approach with expressing the solution as a
% combination of single and double layer potential. 
% The scatterer is a torus with parameters c and a. 

m = 10;
n =10;
c = 2; 
a = 1;
k = 1;
gma = 1;
% The torus is mapped using sheprical harmonics parametrisation to
% rectangle (u, v) \in (0, 2*pi) x (0, 2*pi)

[x1,y1]=ndgrid(2*pi*(1:(m-2))/(m-1),2*pi*(1:n-2)/(n-1)); 
nodes =[x1(:),y1(:)]; % N by 2 matrix listing x,y coordinates of all N=mn nodespi
h1 = 2*pi/ (m-1);
h2 = 2*pi/ (n-1);
N = (m-2)*(n-2);
%position = zeros(N, 1) + 5;
% 1/2 * phi(r)

% Indicate which nodes are boundary, corners. Corners are labelled 0,
% boundaries are labelled 1, 2, 3, 4 anticlockwise with bottom boundary 1.
% Interior nodes are labelled 5.
% for ii = 1:N
%    if (ii==1) | (ii==m) | (ii == N) | (ii == (N -m + 1))
%        position(ii) = 0;
%    elseif mod(ii, m) == 1 
%        position(ii) = 4;
%    elseif mod(ii, m) == 0
%        position(ii) = 2;
%    elseif ii>=1 & ii <= m
%        position(ii) = 1;
%    elseif ii >= (N-m+1) & ii <= N
%        position(ii) = 3
%    else
%        position(ii) = 5;
%           
%    end
%            
% end
% Now we assemble the matrix A and b. b has the boundary values on nodes
% which are known. A has to be constructed using an appropriate quadrature
% scheme. 

b = boundary_val(nodes);
%sigma = gmres(@helmholtz_matvec, b)
 xx = eye(N);
 A = zeros(N, N);
for kk=1:N
    
   A(:, kk) = helmholtz_matvec(xx(:, kk)) ;
        
end

cond(A)
sigma = A\b;
u1 = nodes(5, 1);
v1 = nodes(5, 2);
n_hat = normal_torus(c, a, u1, v1);
M = 1000;
sol = zeros(M, 1);
for eps = (1:M)
r00 = chart(nodes(5, :)) + eps*n_hat;
%r11 = chart(nodes(5, :)) - eps*n_hat;
sol = nonsingular_integrate(r00, sigma, nodes, 0) %+ singular_integrate(r00, x, nodes, 1);
%sol1 = nonsingular_integrate(r11, sigma, nodes, 1) + singular_integrate(r11, x, nodes, 1);
end



    function Ax = helmholtz_matvec(x)
        Ax = zeros(N, 1);
        for ii = 1:N
            r0 = chart(nodes(ii, :));
            si = singular_integrate(ii, x, nodes, 1);
            nsi = nonsingular_integrate(r0, x, nodes, 1);
           Ax(ii) = 0.5*x(ii)  + nsi +si;
                       
        end
    

    end




% for ii=1:N
%    b(ii) = boundary_val(nodes(ii, :));
%    A(ii, ii) = A(ii, ii) + 0.5;
%    A(ii, :) = A(ii, :) + singular_integration(ii, nodes, position);
%    for jj = 1:N
%       A(ii, jj) = A(ii, jj) + nonsingular_kernel(ii, jj, nodes) * ns_weight(ii, jj, position)*h1*h2;
%       %A(ii, jj) = A(ii, jj) + singular_kernel(ii, jj) * s_weight(ii, jj);
%    end    
% end






function r = chart(node)

    u = node(1);
    v = node(2);
    r = [(c + a*cos(v))*cos(u); (c + a*cos(v))*sin(u); a*sin(v)]; 

end

    function nhat = normal_torus(c, a, u, v)
        
        rout = c;
        rin = a;
        iangle = v;
        jangle = u;
        vx = cos(jangle)*(rout+cos(iangle)*rin);
        vy = sin(jangle)*(rout+cos(iangle)*rin);
        vz = sin(iangle)*rin;
        
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
function val_pou = fpou(r0, r1, flag)
    d = 1e-1;    %Check if d is correct
    t = norm(r0 - r1)/ d;
    if t >=1 
        val_pou =0;
    else
    val_pou = exp((2*exp(-1/t))/ (t-1));
    end
    if flag==0
        %disp('hello');
        val_pou =0;
        
    end
end

function val = G1(r0, r1)
    
    val = 1/(4*pi) * cos(k*norm(r0-r1))/ norm(r0-r1);

end


function val = G2(r0, r1)

    val = 1i/(4*pi) * sin(k*norm(r0-r1))/ norm(r0-r1);


end

function val = G3(r0, r1, u, v)
    
    val = 1i/(4*pi) * (sin(k*norm(r0-r1))/ norm(r0-r1)^3  - k*cos(k*norm(r0-r1))/ norm(r0-r1)^2) * sum((r0-r1).* normal_torus(c, a, u, v));


end

function val = G4(r0, r1, u, v)

    val = 1/(4*pi) * (cos(k*norm(r0-r1))/ norm(r0-r1)^3) * sum((r0-r1).* normal_torus(c, a, u, v));


end
function val = G5(r0, r1, u, v)

    val = 1/(4*pi) * (k*sin(k*norm(r0-r1))/ norm(r0-r1)^2) * sum((r0-r1).* normal_torus(c,a, u, v));


end

%Function to calculate the green's function
function val = phi(r0, r1)
    g1 = G1(r0, r1);
    g2 = G2(r0, r1);
    val = g1 + g2;

end

%Function to calculate the normal derivative of green's function 
function val = dphi(r0, r1, u, v)
    g3 = G3 (r0, r1, u, v);
    g4 = G4(r0, r1, u, v);
    g5 = G5(r0, r1, u, v);
    val = g3 + g4 + g5;
end

function val_jac = jacobian(node)
    u = node(1);
    v = node(2);
    val_jac = abs(a*(c + a*cos(v)));

end



%Function for fast interpolation from coarse grid to finer grid
%Returns a matrix that supposedly multiplies the unknown density values at coarse grid and 0 at
%nodes of finer grid to gives values of nodes at finer grid as a linear
%combination of values of  density at coarser grid
% function I = fast_interp(len)
%     %N = m*n;
%     %m_fine = 16*m -15;   %m for finer grid, can change
%     %n_fine = 16*n-15;    %n for finer grid, can change
%     %N_fine = m_fine*n_fine;
%     F = dftmtx(len);
%     Fi = conj(F)/len;
%     [n,m]=size(F);
%     B=permute(reshape(F',m,1,[]),[2 1 3]);
%     p=size(B,3)
%     C=[B;zeros(15,m,p)]
%     F_new=reshape( C(:,:),[],m,1);
%     [n_new, m_new] = size(F_new);
%     F_new = F_new(n_new -15, :);
%     [n,m]=size(F_new');
%     B=permute(reshape(F_new,m,1,[]),[2 1 3]);
%     p=size(B,3)
%     C=[B;zeros(15,m,p)]
%     F_new_T=reshape( C(:,:),[],m,1);
%     F_new = F_new_T';
%     [n_new, m_new] = size(F_new);
%     F_new = F_new(n_new, m_new -15);
%     [m, n] = size(F_new);
%     F_new_actual = dftmtx(m);
%     F_newi = conj(F_new_actual)/m;
%     I = F_newi*F_new;
%     
%     
%     
%     
% 
% end





% %Function to compute nonsingular kernel
% function val = nonsingular_kernel(ii, jj, nodes)
%     source_node = nodes(jj, :);
%     target_node = nodes(ii, :);
%     r1 = chart(source_node); %r1 is the vector of (x, y, z) coordinates of point on surface corresponding to node in parametrised space
%     r0 = chart(target_node);
%     gamma = 1;
%     val = -1i*gamma*(1 - fpou(r0, r1))* phi(r0, r1) * jacobian(source_node)+  (1 - fpou(r0, r1))* dphi(r0, r1) * jacobian(source_node);
% 
% 
% end

function val = nonsingular_integrate( r0, x, nodes, flag)
    %disp('In non singular integrate with '), x 
    val = 0;
    %N = size(nodes, 1);
    %target_node = nodes(ii, :);
    %r0 = chart(target_node);
    for jj= 1:N
        %disp('Iteration over source node '), jj
       source_node = nodes(jj, :);
       if chart(source_node) == r0
           continue
       end
       u = source_node(1);
       v = source_node(2);
       r1 = chart(source_node);
       val = val + h1*h2*(-1i*gma*(1 - fpou(r0, r1, flag))*phi(r0, r1) * jacobian(source_node) * x(jj) ) + h1*h2*(1 - fpou(r0, r1, flag))* dphi(r0, r1, u, v) * jacobian(source_node)*x(jj);
        
    end

end

function node = sing_chart(source_node, target_node)
    rho = source_node(1);
    theta = source_node(2);
    u = target_node(1) + rho*cos(theta);
    v = target_node(2) + rho*sin(theta);
    node = [u, v];

end

	

function Y=interpft2(X,p,q)


X = interpft(X,p);
Y = interpft(X',q);
Y = Y';

end

function xval = get_x(source_node, x)
    u = source_node(1);
    v = source_node(2);
    X = reshape(x, [m-2, n-2])';
    uf = 16;
    X = interpft2(X, uf*(m-2), uf*(n-2));
    mf = uf*(m-2);
    nf = uf*(n-2);
    [xf, yf] = meshgrid(2*pi*(1:mf)/ (mf+1), 2*pi*(1:nf)/ (nf+1));
    xval = interp2(xf, yf, X, u, v, 'spline');
    
    
    
    
    

end


function val = singular_integrate(ii, x, nodes, flag)
    %%disp('in singular integrate with '), x,ii
    val = 0;
    hf1 = h1^2/5;
    hf2 = h2^2/5;
    [x2, y2] = ndgrid(-h1+hf1: hf1:h1-hf1, 0+hf2:hf2:pi-hf2);
    nodes_sing  = [x2(:), y2(:)];
    N_sing = size(nodes_sing, 1);
    
    N = size(nodes, 1);
    target_node = nodes(ii, :);
    r0 = chart(target_node);
    for jj =1 :N_sing
       source_node_sing = nodes_sing(jj, :);
       rho = source_node_sing(1);
       theta = source_node_sing(2);
       source_node = sing_chart(source_node_sing, target_node);
       u = source_node(1);
       v = source_node(2);
       r1 = chart(source_node);
       xval = get_x(source_node, x);
       val = val +  hf1*hf2*(-1i*gma*(fpou(r0, r1, flag))*phi(r0, r1) * jacobian(source_node) * xval * abs(rho) ) ;
       val = val + hf1*hf2*(fpou(r0, r1, flag))* dphi(r0, r1, u, v) * jacobian(source_node)*xval*abs(rho);
        
    end

end


% %Function to compute weight associated to node in 2D trapezoidal rule
% function weight = ns_weight(ii, jj, position)
%     if position(jj) == 0  % corner 
%         weight = 0.25;
%     elseif position(jj) == 1 | position(jj) == 2 | position(jj) == 3 | position(jj) == 4  % boundary
%         weight = 0.5; 
%         
%     else                    % Interior node
%         weight = 1;    
%     end
%     
% end


% Function that gives the given boundary Dirchlet data
function val = boundary_val(nodes)
    u = nodes(:, 1);
    v = nodes(:, 2);
    val = ones(size(u, 1), 1) ;
end


end