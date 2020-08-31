function projectors1Qx = Projectors(N)
    %UNTITLED Summary of this function goes here
    %   Detailed explanation goes here
    %Generate 6^n Pauli kets resulting from all kronecker products
    %of the 6 Pauli's eigenvectors 
    
    h = [1 ; 0];
    v = [0 ; 1];
    d = (h + v)/sqrt(2);
    a = (h - v)/sqrt(2);
    r = (h + v*1i)/sqrt(2);
    l = (h - v*1i)/sqrt(2);
    
    projectors1Q = zeros(2,2,6);
    projectors1Q(:,:,1) = h*h';
    projectors1Q(:,:,2) = v*v';
    projectors1Q(:,:,3) = d*d';
    projectors1Q(:,:,4) = a*a';
    projectors1Q(:,:,5) = r*r';
    projectors1Q(:,:,6) = l*l';
    projectors1Qx = projectors1Q;
    
    
    projectorsNQs = zeros(2^N,2^N,6^N);
    projNQs = zeros(4,4,36);
    
    if N > 1
        for i=2:N
           projNQs = zeros(2^(i),2^(i),6^(i));
           c = 0;
           size1 = size(projectors1Qx);
           size1 = size1(3);
           size2 = size(projectors1Q);
           size2 = size2(3);
           for j=1:size1
               for k=1:size2
                c = c + 1;
                projNQs(:,:,c) = kron(projectors1Qx(:,:,j), projectors1Q(:,:,k));
               end
           end
           projectors1Qx = projNQs;
           %projNQs = zeros(2^(i + 1),2^(i + 1),6^(i + 1));
        end
    
end
