close all; 
clear; 
clc; 


xs = im2double(imread('building_roof.jpg'));

f = @(x) imbilatfilt(x,0.05,3);

ys = f(xs);

maxiter = 10;

Tc = ys;
Sc = ys;
Pc = ys;

for c=1:size(xs,3)
    yc = ys(:,:,c);
    
    T = yc;
    S = yc;
    P = yc;
    
    for n=1:maxiter        
        ht = yc-f(T);
        hs = yc-f(S);
        hp = yc-f(P);
        
        % S-method
        d = f(S+hs)-f(S);
        S = S + norm(hs)*hs/norm(d);
        
        % P-method
        d = (f(P+hp)-f(P-hp))/2;
        lam = norm(hp)^2/(norm(d)+eps)^2;
        P = P + lam*d;
        
        % Tao et al.
        T = T + ht;
        
    end

    Sc(:,:,c) = S;
    Pc(:,:,c) = P;
    Tc(:,:,c) = T;
end

figure,imshow([xs,ys]),title('Initial and filtered image');
figure,imshow([Tc,Sc,Pc]),title('T  S  P');

