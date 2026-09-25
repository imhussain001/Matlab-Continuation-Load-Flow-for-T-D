function sys = ieee33_data()
%IEEE33_DATA  IEEE 33-bus radial distribution feeder (Baran & Wu, 1989).
%
%   sys = ieee33_data() returns
%     sys.fromBus, sys.toBus  32 x 1 branch terminals (bus 1 = substation)
%     sys.z                   32 x 1 complex branch impedance (pu)
%     sys.Sload               33 x 1 complex base load (pu), load positive
%     sys.BusData/BranchData  same system in the IEEE-14 matrix format so
%                             that run_cpf (Newton CPF) can be applied
%     sys.Vbase_kV = 12.66, sys.Sbase_MVA = 100, sys.Zbase_ohm
%
%   Raw data: R, X in ohm; P in MW, Q in Mvar. Total load 3.715 MW +
%   j2.300 Mvar.
%
%   NOTE: bus 33 reactive load restored to the standard 0.040 Mvar. The
%   legacy main.m had -10.100 Mvar - a 10 Mvar capacitor on a 3.7 MW
%   feeder - which is not Baran & Wu data.

Vbase_kV = 12.66;
Sbase = 100;
Zbase = Vbase_kV^2 / Sbase;

%        id from to   R(ohm)  X(ohm)
branch = [1   1   2   0.0922  0.0470
          2   2   3   0.4930  0.2511
          3   3   4   0.3660  0.1864
          4   4   5   0.3811  0.1941
          5   5   6   0.8190  0.7070
          6   6   7   0.1872  0.6188
          7   7   8   0.7114  0.2351
          8   8   9   1.0300  0.7400
          9   9  10   1.0040  0.7400
         10  10  11   0.1996  0.0650
         11  11  12   0.3744  0.1238
         12  12  13   1.4680  1.1550
         13  13  14   0.5416  0.7129
         14  14  15   0.5910  0.5260
         15  15  16   0.7463  0.5450
         16  16  17   1.2890  1.7210
         17  17  18   0.7320  0.5740
         18   2  19   0.1640  0.1565
         19  19  20   1.5042  1.3554
         20  20  21   0.4095  0.4784
         21  21  22   0.7089  0.9373
         22   3  23   0.4512  0.3083
         23  23  24   0.8980  0.7091
         24  24  25   0.8960  0.7011
         25   6  26   0.2030  0.1034
         26  26  27   0.2842  0.1447
         27  27  28   1.0590  0.9337
         28  28  29   0.8042  0.7006
         29  29  30   0.5075  0.2585
         30  30  31   0.9744  0.9630
         31  31  32   0.3105  0.3619
         32  32  33   0.3410  0.5302];

%     bus  P(MW)  Q(Mvar)
loadData = [1  0.000  0.000
        2  0.100  0.060
        3  0.090  0.040
        4  0.120  0.080
        5  0.060  0.030
        6  0.060  0.020
        7  0.200  0.100
        8  0.200  0.100
        9  0.060  0.020
       10  0.060  0.020
       11  0.045  0.030
       12  0.060  0.035
       13  0.060  0.035
       14  0.120  0.080
       15  0.060  0.010
       16  0.060  0.020
       17  0.060  0.020
       18  0.090  0.040
       19  0.090  0.040
       20  0.090  0.040
       21  0.090  0.040
       22  0.090  0.040
       23  0.090  0.050
       24  0.420  0.200
       25  0.420  0.200
       26  0.060  0.025
       27  0.060  0.025
       28  0.060  0.020
       29  0.120  0.070
       30  0.200  0.600
       31  0.150  0.070
       32  0.210  0.100
       33  0.060  0.040];

N = size(loadData,1);
L = size(branch,1);

sys.name = 'IEEE 33-bus';
sys.Vbase_kV = Vbase_kV;
sys.Sbase_MVA = Sbase;
sys.Zbase_ohm = Zbase;
sys.fromBus = branch(:,2);
sys.toBus   = branch(:,3);
sys.z       = complex(branch(:,4), branch(:,5)) / Zbase;
sys.Sload   = complex(loadData(:,2), loadData(:,3)) / Sbase;

% Same system in BusData/BranchData format (see data/ieee14_data.m)
BusData = zeros(N,17);
BusData(:,1) = (1:N).';
BusData(:,2:3) = 1;
BusData(:,4) = 0;  BusData(1,4) = 3;          % bus 1 = slack
BusData(:,5) = 1;                              % slack |V| = 1.0 pu
BusData(:,7) = real(sys.Sload);
BusData(:,8) = imag(sys.Sload);
BranchData = zeros(L,21);
BranchData(:,1) = sys.fromBus;
BranchData(:,2) = sys.toBus;
BranchData(:,3:5) = 1;
BranchData(:,7) = real(sys.z);
BranchData(:,8) = imag(sys.z);
sys.BusData = BusData;
sys.BranchData = BranchData;
end
