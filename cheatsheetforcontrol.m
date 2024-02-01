clc, clear, close

m = 1;
k = 1;
b = 0.2;
F = 1;

A = [0 1; -k/m -b/m];
B = [0 1/m]';
C = [1 0];
D = [0];

sys = ss(A,B,C,D)

% Define system by TF manual input
s = tf('s');
sys = 1/(m*s^2+b*s+k)

% Define system by TF num, den
num = [1];
den = [m b k];
sys = tf(num,den)

% Get SS from TF
s = tf('s');
G = (2*s+1)/(4*s^2+3*s+2)
[A,B,C,D] = ssdata(G)

% Get num,den from SS
[num,den] = ss2tf(A,B,C,D)

% Get SS from num,den
[A,B,C,D] = tf2ss(num,den)

% Get zeros poles k from num,den 
num=[0 0.5 0.25]
den=[1 0.75 0.5]
[z,p,k] = tf2zp(num,den)

% Get zeros poles k from SS 
A =[-0.7500   -0.5000;1.0000 0]   
B =[1;0]
C=[0.5000    0.2500]
D =[0]
[z,p,k] = ss2zp(A,B,C,D)

% Get SS from zeros poles k
[A,B,C,D] = zp2ss(z,p,k)

% Get num,den from zeros poles k
[num,den] = zp2tf(z,p,k)

% Get Partial Fraction Expansion(PFE) of the given TF
num = [1 2 1]
den = [1 2 3 4 5]
[r,p,k] = residue(num, den)
% k is non empty only when OrdNum > OrdDen

% Controllability and Observability
num = [1 12];
den = [1 4 4 2];
% sys = tf(num,den)

[A, B, C, D] = tf2ss(num,den);
Co = ctrb(A,B);
rank(Co);
% uncontrollable part
unco = length(A) - rank(Co)

Ob = obsv(A,C);
rank(Ob);
%unobservable part
unobsv = length(A) - rank(Ob)

% Polynomial related funcs
P1 = [3 0 0 4 77 8];
P2 = [5 4 2];
polyval(P1,3)
P3 = conv(P1,P2)
deconv(P3,P2)

% Matrix related funcs
B = [1 9 4; 2 3 4; 1 2 7];
diag(B)
trace(B)
[V,D] = eig(B);
[U,S,D] = svd(B);
eye(3)
ones(3)
zeros(5)

% others
vpa(pi,15)
floor(6.2)
