clear
clc

%symbolic test for question 1 and 2
syms sl1 sl2 sl3 st1 st2 st3 std1 std2 std3 stdd1 stdd2 stdd3 sm1 sm2 sm3 sfx sfy sfz stx sty stz sml

%Uncomment to view full symbolic derivation 
%hw4(sl1, sl2, sl3, st1, st2, st3, std1, std2, std3, stdd1, stdd2, stdd3, sm1, sm2, sm3, sfx, sfy, sfz, stx, sty, stz, sml, 1);

%numerical test
m1 = 2;
m2 = 2;
m3 = 2;
ml = 1;

l1 = 50; 
l2 = 25;
l3 = 10;

t1 = 35;
t2 = 20;
t3 = 30;

td1 = 45;
td2 = 45;
td3 = 45;

tdd1 = 0;
tdd2 = 0;
tdd3 = 0;

fx = 0;
fy = ml*-9.8;
fz = 0;

tx = 0;
ty = 0;
tz = 0;

%Uncomment for question 3 and 4
%hw4(l1, l2, l3, t1, t2, t3, td1, td2, td3, tdd1, tdd2, tdd3, m1, m2, m3, fx, fy, fz, tx, ty, tz, ml, 0); 
hw4(l1, l2, l3, st1, st2, st3, td1, td2, td3, tdd1, tdd2, tdd3, m1, m2, m3, fx, fy, fz, tx, ty, tz, ml, 0); 