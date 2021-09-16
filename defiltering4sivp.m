close all;
clear all;
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
    xc = xs(:,:,c);
    
    nys = norm(yc(:));
    nxs = norm(xc(:));
    
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
        
        te = norm(xc(:)-T(:))/nxs;
        se = norm(xc(:)-S(:))/nxs;
        pe = norm(xc(:)-P(:))/nxs;
        
        fprintf('%d %d %f %f %f\n',[c,n,te,se,pe]);
        
    end
    
    Tc(:,:,c) = T;
    Pc(:,:,c) = P;
    Sc(:,:,c) = S;
    
end

figure,imshow([xs,ys]);
figure,imshow([Tc,Sc,Pc]),title('T  S  P');
