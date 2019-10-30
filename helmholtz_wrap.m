function [] = helmholtz_wrap()
% k = 3;
% for m = 8:8:32
%     helmholtz_sphereVec(m, m, k, 1);
% end

k = 0;

for m = 4:4:32
    helmholtz_sphereVec(m, m, k, 0);
end

end