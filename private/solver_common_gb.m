% The common Groebner basis solver for the 1st sub-problem
function R = solver_common_gb(M)

    sols = solver_vpnp(M);
    sols = sols(:, any(not(imag(sols))));
    
    
    r11 = sols(1, :);
    r12 = sols(2, :);
    r13 = sols(3, :);
    r21 = sols(4, :);
    r22 = sols(5, :);
    r23 = 1;

    nsols = length(r11);
    R1    = zeros(3,3,nsols);
    R2    = zeros(3,3,nsols);
    for i=1:nsols
        r1 = [r11(i), r12(i), r13(i)];
        r1 = r1 / norm(r1);
        r2 = [r21(i), r22(i), r23];
        r2 = r2 / norm(r2);
        r3 = cross(r1, r2);

        R_tmp = [r1
                 r2
                 r3];
                 
        % (optional) Enforce R as a rotation matrix
        % Not necessary if root polishing is enabled.
        %[U,~,V] = svd( R_tmp );
        %R_tmp   = U * diag([1,1,det(U*V')]) * V';
        
        R1(:,:,i) = R_tmp;
        R2(:,:,i) = diag([-1,-1,1]) * R_tmp; % == [-r1; -r2; r3]
    end
    R = cat(3,R1,R2);

end


% The following code is generated by using Larsson's automatic generator
% for polynomial solvers v0.5. 
% To reproduce the code, install the generator and run
%   opts = default_options();
%   solv = generate_solver('vpnp', @problem_vpnp, opts);
function sols = solver_vpnp(data)
[C0,C1] = setup_elimination_template(data);
C1 = C0 \ C1;
RR = [-C1(end-8:end,:);eye(20)];
AM_ind = [27,12,13,1,16,2,17,3,21,4,5,22,6,25,7,26,8,28,29,9];
AM = RR(AM_ind,:);
[V,D] = eig(AM);
V = V ./ (ones(size(V,1),1)*V(1,:));
sols(1,:) = V(2,:);
sols(2,:) = V(5,:);
sols(3,:) = V(9,:);
sols(4,:) = V(14,:);
sols(5,:) = diag(D).';
end

% Action =  x5
% Quotient ring basis (V) = 1,x1,x1*x5,x1*x5^2,x2,x2*x4,x2*x5,x2*x5^2,x3,x3^2,x3*x4,x3*x5,x3*x5^2,x4,x4^2,x4*x5,x4*x5^2,x5,x5^2,x5^3,
% Available monomials (RR*V) = x1*x5^3,x2*x4*x5,x2*x5^3,x3^2*x5,x3*x4*x5,x3*x5^3,x4^2*x5,x4*x5^3,x5^4,1,x1,x1*x5,x1*x5^2,x2,x2*x4,x2*x5,x2*x5^2,x3,x3^2,x3*x4,x3*x5,x3*x5^2,x4,x4^2,x4*x5,x4*x5^2,x5,x5^2,x5^3,
function [coeffs] = compute_coeffs(data)
coeffs(1) = data(3);
coeffs(2) = data(9);
coeffs(3) = -data(2);
coeffs(4) = -data(8) + data(15);
coeffs(5) = -data(14);
coeffs(6) = data(21);
coeffs(7) = -data(20);
coeffs(8) = data(31);
coeffs(9) = data(27) + data(32);
coeffs(10) = -data(26) + data(33);
coeffs(11) = data(24);
coeffs(12) = data(30);
coeffs(13) = -data(25);
coeffs(14) = -data(27) - data(32);
coeffs(15) = -data(23);
coeffs(16) = -data(29) + data(36);
coeffs(17) = -data(35);
coeffs(18) = -data(3);
coeffs(19) = -data(9);
coeffs(20) = data(1) - data(15);
coeffs(21) = data(7);
coeffs(22) = data(13);
coeffs(23) = -data(21) - data(31);
coeffs(24) = -data(32);
coeffs(25) = data(19) - data(33);
coeffs(26) = -data(24);
coeffs(27) = -data(27);
coeffs(28) = data(25);
coeffs(29) = -data(30);
coeffs(30) = data(20);
coeffs(31) = data(21) + data(31);
coeffs(32) = data(22) - data(36);
coeffs(33) = data(28);
coeffs(34) = data(34);
coeffs(35) = data(2);
coeffs(36) = -data(1) + data(8);
coeffs(37) = -data(7);
coeffs(38) = data(14);
coeffs(39) = -data(13);
coeffs(40) = data(20) + data(25);
coeffs(41) = -data(19) + data(26);
coeffs(42) = data(27);
coeffs(43) = data(23);
coeffs(44) = -data(20) - data(25);
coeffs(45) = -data(21);
coeffs(46) = -data(22) + data(29);
coeffs(47) = -data(28);
coeffs(48) = data(32);
coeffs(49) = -data(31);
coeffs(50) = data(35);
coeffs(51) = -data(34);
coeffs(52) = -data(19);
coeffs(53) = -data(26);
coeffs(54) = -data(33);
coeffs(55) = data(1) - data(22);
coeffs(56) = data(7) - data(23);
coeffs(57) = data(13) - data(24);
coeffs(58) = data(19);
coeffs(59) = data(2) - data(28);
coeffs(60) = data(8) - data(29);
coeffs(61) = data(14) - data(30);
coeffs(62) = data(26);
coeffs(63) = data(3) - data(34);
coeffs(64) = data(9) - data(35);
coeffs(65) = data(15) - data(36);
coeffs(66) = data(33);
coeffs(67) = 1;
coeffs(68) = -1;
end
function [C0,C1] = setup_elimination_template(data)
[coeffs] = compute_coeffs(data);
coeffs0_ind = [18,35,52,67,19,1,36,44,2,37,53,67,18,35,52,67,20,3,38,19,1,36,44,23,21,4,39,2,37,53,14,67,20,3,38,23,22,5,21,4,39,14,54,67,22,5,54,67,18,35,...
52,67,19,1,35,18,36,52,44,67,2,1,36,19,37,44,53,67,2,37,53,67,20,3,18,38,52,23,35,67,23,40,55,21,4,3,38,19,20,39,44,23,14,36,1,24,6,41,56,4,...
39,21,53,14,37,67,2,23,40,55,22,5,20,23,54,38,3,67,25,7,42,24,6,41,56,57,5,21,22,14,54,39,4,67,25,7,42,57,22,54,67,5,23,40,55,18,35,52,67,67,...
24,6,40,23,41,55,56,19,36,1,44,67,6,41,24,56,37,2,53,67,25,7,23,42,55,57,20,38,40,3,23,67,26,43,58,7,42,24,25,56,57,21,39,41,4,14,6,68,26,43,...
58,68,25,57,22,42,5,54,67,7,26,43,58,23,40,55,68,67,43,26,58,24,41,6,56,68,26,58,25,42,43,7,57,68,26,43,58,68,18,35,52,67,19,35,36,1,52,44,67,18,...
36,1,37,2,44,53,19,67,37,2,53,67,20,38,35,3,52,23,18,67,27,8,41,59,21,38,3,39,36,4,44,23,1,14,19,20,9,44,60,39,4,37,53,14,2,67,21,27,8,41,...
59,22,38,5,23,3,54,20,67,28,10,45,9,44,60,61,5,39,14,54,4,21,67,22,28,10,45,61,54,5,22,67,27,8,41,59,23,40,55,18,52,35,67,67,9,8,41,27,44,59,...
60,24,40,41,6,55,56,19,1,44,36,23,67,67,9,44,60,41,6,56,2,53,37,67,24,67,28,10,27,45,59,61,41,8,25,42,40,7,55,57,20,3,23,38,23,67,29,11,46,40,...
10,45,28,60,61,44,9,42,7,41,56,57,6,21,4,14,39,24,25,67,29,11,46,40,28,61,45,10,42,57,7,22,5,54,67,25,29,11,46,40,27,41,8,59,26,43,58,23,55,40,...
68,67,11,46,29,40,44,9,60,43,58,24,6,56,41,68,26,67,29,40,28,45,46,10,61,11,43,58,25,7,57,42,26,68,29,46,11,40,26,58,43,68,27,41,8,59,52,35,67,18,...
41,8,44,9,59,60,27,44,36,19,1,67,44,9,60,53,37,67,2,67,28,45,41,10,59,8,61,27,23,38,20,3,12,47,62,68,45,10,44,60,61,9,28,14,39,21,4,67,12,47,...
62,68,45,61,10,28,54,67,22,5,12,47,62,29,46,11,40,27,8,59,41,55,40,23,68,67,12,47,62,68,46,11,40,9,60,44,29,56,41,24,6,67,62,47,68,12,46,40,11,28,...
10,61,45,29,57,42,25,7,47,12,62,68,29,11,40,46,58,43,68,26,52,18,67,35,44,52,35,67,19,18,1,36,53,44,1,36,19,2,67,37,53,2,37,67,23,67,20,3,38,35,...
18,52,25,13,48,63,14,23,3,38,1,21,20,4,39,36,19,44,30,10,49,64,14,4,39,67,2,21,37,53,25,13,48,63,54,3,22,5,67,38,20,23,31,14,30,10,49,64,65,54,...
5,4,67,22,39,21,14,31,14,65,67,5,22,54,25,13,48,63,55,23,67,35,18,40,52,67,30,10,13,48,25,49,63,64,56,55,40,24,23,6,1,36,19,41,44,67,10,49,30,64,...
56,6,41,24,67,2,37,53,31,14,25,63,65,48,13,57,25,7,3,38,20,42,23,40,23,55,67,67,32,15,50,31,14,30,31,64,65,49,10,57,7,42,6,25,4,39,21,14,41,24,...
56,67,32,15,50,31,31,65,14,7,67,5,22,54,42,25,57,67,32,15,50,31,25,48,13,63,58,26,40,23,68,43,55,67,15,50,32,31,30,49,10,64,58,43,68,26,6,41,24,56,...
32,31,31,50,14,65,15,67,68,7,42,25,57,43,26,58,32,50,15,31,68,43,26,58,25,48,13,63,59,27,8,41,52,18,67,35,30,48,13,49,10,63,64,25,60,59,8,41,27,9,...
44,44,19,1,36,67,49,10,64,30,60,9,44,53,2,67,37,67,31,48,14,63,13,65,25,61,8,28,10,45,41,27,23,20,3,59,38,67,33,16,51,9,14,49,64,65,10,30,31,61,...
10,45,9,28,44,14,21,4,60,39,67,67,33,16,51,9,32,50,15,31,25,13,63,48,40,29,11,8,41,27,46,59,55,23,67,40,63,25,13,48,67,52,35,18,64,63,13,48,30,25,...
10,49,1,44,36,19,64,10,49,30,67,2,53,37,65,13,31,14,48,25,63,3,23,38,20,67,34,17,66,68,65,14,10,31,49,30,64,4,14,39,21,67,34,17,66,68,31,32,15,13,...
48,25,50,63,55,40,23,67,47,12,62,68,59,41,27,8,16,51,33,9,50,15,31,30,10,64,49,32,40,11,46,29,9,44,60,67,56,24,6,41,47,12,62,68,60,44,9,67,33,16,...
51,9,65,14,31,67,10,45,28,54,22,5,67,61,33,9,51,16,50,31,15,31,14,65,32,67,11,10,45,28,61,46,29,57,25,7,40,42,47,62,12,68,61,45,28,10,33,51,16,9,...
32,15,31,50,11,46,29,40,58,26,68,43,12,62,47,68,40,46,29,11,62,47,68,12];
coeffs1_ind = [68,17,66,34,66,34,17,68,13,63,48,25,34,17,66,68,9,33,16,51,63,25,13,8,59,41,27,48,33,51,16,9,63,48,25,13,62,12,68,47,59,27,8,41,66,17,68,34,10,64,...
49,30,17,34,66,68,31,15,50,32,10,49,30,64,6,56,41,24,17,66,68,34,9,16,51,33,64,30,10,9,60,44,67,49,51,16,9,33,64,49,30,10,62,12,47,68,60,9,67,44,...
68,17,34,66,14,65,31,67,34,17,66,68,14,31,65,67,5,54,22,67,34,66,68,17,15,14,31,65,50,67,32,31,7,57,42,25,66,17,34,68,16,51,33,65,31,14,9,67,10,61,...
45,28,51,9,16,33,65,31,14,67,68,12,47,61,28,10,62,45,68,17,34,66,15,31,50,32,34,17,66,68,15,50,32,31,68,58,43,26,34,17,66,68,16,51,33,9,31,32,15,11,...
40,46,29,50,33,16,9,51,31,50,32,15,68,12,47,62,40,29,11,46,66,34,17,68,16,9,51,33,66,68,34,17,9,33,16,68,12,62,47,51,9,51,33,16,62,12,68,47];
C0_ind = [1,3,8,29,101,102,103,108,202,203,208,229,304,306,307,309,401,402,403,404,405,406,407,408,501,502,503,505,506,507,508,509,604,605,606,607,701,702,704,705,706,707,708,729,804,805,807,809,910,916,...
919,955,1010,1011,1013,1015,1016,1018,1019,1028,1111,1112,1113,1115,1116,1118,1119,1155,1212,1213,1218,1228,1310,1311,1314,1316,1317,1319,1322,1326,1401,1403,1408,1410,1411,1412,1413,1414,1415,1416,1417,1418,1419,1422,1427,1501,1502,1503,1508,1512,...
1513,1515,1517,1518,1522,1526,1527,1604,1606,1607,1610,1611,1614,1617,1619,1622,1627,1655,1701,1702,1703,1704,1705,1706,1707,1708,1712,1714,1715,1717,1718,1722,1727,1728,1804,1805,1806,1807,1814,1817,1826,1827,1910,1916,1919,1920,1921,1924,1925,2000,...
2010,2011,2013,2015,2016,2018,2019,2020,2021,2023,2024,2097,2112,2113,2115,2118,2121,2123,2124,2125,2210,2211,2214,2216,2217,2219,2220,2221,2222,2223,2224,2292,2301,2303,2308,2312,2313,2314,2315,2317,2318,2320,2321,2322,2323,2324,2327,2329,2404,2406,...
2407,2409,2414,2417,2420,2422,2423,2424,2425,2427,2510,2516,2519,2520,2521,2524,2555,2558,2613,2615,2618,2620,2621,2623,2624,2628,2714,2717,2720,2721,2722,2723,2724,2726,2820,2821,2824,2825,2930,2933,2939,2949,3030,3031,3033,3035,3037,3039,3047,3048,...
3131,3132,3133,3135,3137,3139,3148,3149,3231,3232,3237,3247,3330,3333,3334,3335,3336,3339,3345,3346,3401,3402,3403,3408,3430,3431,3432,3433,3434,3435,3436,3437,3438,3439,3445,3448,3502,3503,3508,3531,3532,3534,3536,3537,3538,3546,3548,3604,3605,3606,...
3607,3630,3634,3635,3636,3638,3639,3645,3649,3701,3702,3703,3705,3706,3707,3708,3732,3734,3736,3737,3738,3745,3747,3748,3804,3805,3806,3807,3836,3838,3845,3846,3910,3911,3916,3919,3930,3933,3939,3940,3942,3943,3944,3999,4011,4012,4013,4015,4016,4018,...
4019,4030,4031,4033,4035,4037,4039,4040,4041,4042,4043,4048,4094,4100,4112,4113,4118,4131,4132,4137,4141,4142,4143,4144,4148,4197,4210,4211,4214,4216,4217,4219,4222,4227,4230,4233,4234,4235,4236,4239,4240,4241,4242,4243,4245,4259,4301,4302,4303,4308,...
4312,4313,4315,4317,4318,4322,4327,4331,4332,4334,4336,4337,4338,4340,4341,4342,4343,4345,4348,4392,4404,4405,4406,4407,4414,4417,4422,4427,4434,4436,4438,4440,4441,4442,4444,4445,4510,4511,4516,4519,4520,4521,4523,4524,4530,4533,4539,4540,4542,4543,...
4549,4557,4612,4613,4615,4618,4621,4623,4624,4631,4637,4640,4641,4642,4643,4647,4648,4658,4714,4717,4720,4721,4722,4723,4724,4727,4734,4736,4740,4741,4742,4743,4745,4746,4820,4821,4823,4824,4840,4842,4843,4844,4930,4933,4935,4939,4950,4951,4952,4953,...
5031,5032,5033,5035,5037,5039,5048,5050,5051,5053,5054,5099,5131,5132,5137,5150,5151,5152,5154,5194,5230,5233,5234,5235,5236,5238,5239,5245,5250,5251,5253,5254,5302,5303,5308,5329,5331,5332,5334,5336,5337,5338,5348,5350,5351,5353,5354,5359,5405,5406,...
5407,5409,5434,5436,5438,5445,5450,5452,5453,5454,5511,5516,5519,5530,5533,5535,5539,5540,5541,5542,5543,5550,5551,5553,5555,5556,5612,5613,5618,5628,5631,5632,5637,5641,5642,5643,5648,5650,5651,5653,5654,5657,5717,5722,5726,5727,5734,5736,5738,5740,...
5741,5742,5743,5745,5750,5751,5753,5754,5821,5823,5824,5825,5840,5841,5842,5843,5850,5851,5852,5853,5960,5967,5974,5975,6060,6061,6063,6066,6067,6068,6069,6075,6160,6161,6162,6163,6168,6169,6174,6175,6261,6262,6263,6266,6360,6364,6367,6369,6375,6377,...
6379,6384,6401,6402,6403,6408,6460,6461,6462,6463,6465,6467,6468,6469,6475,6477,6479,6484,6501,6502,6503,6508,6561,6562,6563,6564,6565,6568,6577,6584,6604,6605,6606,6607,6660,6665,6667,6669,6674,6677,6679,6684,6701,6702,6704,6705,6706,6707,6708,6761,...
6762,6765,6766,6768,6777,6779,6784,6804,6805,6807,6864,6865,6879,6884,6910,6911,6916,6919,6960,6967,6970,6972,6973,6975,6976,6998,7010,7011,7012,7013,7015,7016,7018,7019,7060,7061,7063,7067,7068,7069,7071,7072,7073,7075,7076,7095,7112,7113,7115,7118,...
7161,7162,7163,7168,7170,7171,7172,7176,7210,7211,7214,7217,7219,7222,7227,7260,7267,7269,7271,7272,7273,7275,7276,7277,7279,7284,7296,7300,7301,7302,7303,7308,7312,7314,7315,7317,7318,7322,7327,7361,7362,7363,7365,7368,7371,7372,7373,7376,7377,7379,...
7384,7397,7404,7405,7406,7407,7414,7417,7427,7465,7470,7471,7473,7476,7477,7479,7484,7492,7510,7511,7516,7519,7520,7521,7523,7524,7560,7567,7572,7573,7574,7575,7576,7578,7612,7613,7615,7618,7620,7621,7623,7624,7661,7663,7666,7668,7671,7672,7673,7676,...
7714,7717,7720,7722,7723,7724,7727,7758,7764,7771,7772,7773,7776,7777,7779,7784,7820,7821,7823,7824,7870,7872,7873,7876,7930,7933,7935,7939,7960,7967,7969,7975,7980,7981,7983,7993,8030,8031,8032,8033,8035,8037,8039,8048,8060,8061,8062,8063,8068,8069,...
8075,8080,8081,8082,8093,8098,8131,8132,8137,8148,8161,8162,8163,8180,8182,8183,8193,8195,8230,8234,8235,8236,8238,8239,8245,8260,8265,8267,8269,8275,8277,8279,8280,8281,8282,8284,8293,8299,8301,8302,8303,8308,8332,8334,8336,8337,8338,8345,8348,8361,...
8362,8363,8365,8368,8377,8380,8381,8382,8384,8393,8394,8396,8410,8411,8416,8419,8430,8433,8435,8439,8440,8441,8442,8443,8460,8467,8469,8471,8472,8473,8475,8476,8480,8481,8485,8493,8560,8567,8569,8575,8586,8588,8589,8590,8660,8661,8662,8663,8667,8668,...
8669,8675,8687,8688,8689,8690,8761,8762,8763,8768,8786,8787,8788,8789,8860,8865,8867,8869,8877,8879,8884,8887,8888,8889,8890,8898,8901,8902,8908,8929,8961,8962,8965,8968,8977,8979,8984,8987,8988,8989,8990,8995,9010,9011,9019,9055,9060,9067,9069,9071,...
9072,9073,9075,9076,9088,9089,9090,9091,9133,9135,9139,9149,9150,9151,9153,9154,9212,9213,9215,9218,9231,9232,9237,9240,9241,9242,9243,9248,9261,9262,9263,9268,9271,9272,9276,9278,9280,9281,9282,9293,9331,9332,9337,9347,9350,9351,9354,9356,9404,9405,...
9406,9407,9436,9438,9445,9459,9465,9477,9479,9480,9481,9482,9483,9484,9514,9517,9522,9527,9534,9536,9538,9540,9541,9542,9545,9557,9565,9571,9572,9573,9576,9577,9579,9580,9581,9582,9584,9593,9634,9636,9638,9646,9650,9651,9653,9654,9720,9721,9723,9724,...
9740,9741,9742,9743,9771,9772,9773,9776,9780,9781,9783,9793,9841,9842,9843,9844,9850,9851,9853,9854,9950,9951,9952,9954];
C1_ind = [86,87,88,90,160,167,169,174,187,188,189,190,230,235,239,249,260,267,269,275,280,281,282,287,288,289,290,293,330,333,335,339,350,351,353,354,360,369,374,375,380,381,382,393,461,462,466,468,487,488,...
489,490,512,515,518,528,561,562,563,568,571,572,573,576,587,588,589,590,632,637,647,648,661,662,663,668,680,681,682,687,688,689,691,693,731,732,737,748,750,751,753,754,761,762,763,766,780,782,785,793,...
864,865,879,884,887,888,890,891,904,905,907,909,965,979,984,986,987,988,990,996,1014,1017,1026,1027,1065,1071,1073,1076,1077,1078,1079,1084,1087,1088,1089,1090,1136,1138,1145,1146,1165,1177,1179,1180,1181,1182,1184,1185,1187,1188,...
1189,1190,1234,1236,1238,1245,1250,1253,1254,1256,1264,1265,1277,1280,1281,1282,1284,1293,1370,1371,1373,1376,1387,1388,1389,1390,1420,1423,1424,1425,1471,1472,1473,1476,1486,1488,1489,1490,1540,1541,1542,1544,1571,1572,1573,1576,1580,1581,1582,1587,...
1588,1589,1590,1593,1640,1641,1642,1643,1650,1651,1653,1654,1670,1671,1672,1676,1680,1681,1682,1693,1780,1781,1782,1783,1787,1788,1789,1790,1850,1852,1853,1854,1880,1881,1882,1886,1887,1888,1889,1893,1950,1951,1953,1954,1980,1982,1983,1993];
C0 = zeros(100,100);
C1 = zeros(100,20);
C0(C0_ind) = coeffs(coeffs0_ind);
C1(C1_ind) = coeffs(coeffs1_ind);
end