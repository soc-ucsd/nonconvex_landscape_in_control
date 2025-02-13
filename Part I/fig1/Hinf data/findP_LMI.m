function [P, tmin] = findP_LMI(A, B, C, D, hinf_norm)
[nx, ~] = size(A);

setlmis([]);
Pvar = lmivar(1, [nx, 1]);
lmiterm([1 1 1 Pvar], A', 1, 's');
lmiterm([1 1 2 Pvar], 1, B);
lmiterm([1 1 3 0], C');
lmiterm([1 2 2 0], -(hinf_norm));
lmiterm([1 2 3 0], D');
lmiterm([1 3 3 0], -(hinf_norm));
lmisys = getlmis;

[tmin, xfeas] = feasp(lmisys, [0 0 -1 0 1]);
P = dec2mat(lmisys, xfeas, Pvar);
end