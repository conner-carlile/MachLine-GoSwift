import pickle
import csv
import numpy as np
from off_body import *
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.animation import PillowWriter
from scipy.signal import savgol_filter
from scipy.optimize import fsolve
from exponentalSmoothing import exponential_smoothing
from rapidboom.sboomwrapper import SboomWrapper
import pyldb

## Read in a CSV file with 2 columns of data x,p
#x_loc = []
with open('studies/Goswift/results/pod_wedge_x_loc.csv', 'r') as csvfile:
    reader = csv.reader(csvfile)
    next(reader)  # Skip the header row if it exists
    x_loc = [float(row[0]) for row in reader]

with open('studies/Goswift/results/pod_wedge_nearfield.csv', 'r') as csvfile:
    reader = csv.reader(csvfile)
    next(reader)  # Skip the header row if it exists
    nearfield_sig = [float(row[0]) for row in reader]
#print("x_loc: ", x_loc)


sav = savgol_filter(nearfield_sig, 100, 3)
exp = exponential_smoothing(nearfield_sig, 0.05)

#plt.plot(x_loc, nearfield_sig)
#plt.title('Raw Data')
#plt.show()
#plt.plot(x_loc, sav)
#plt.title('Savitzky-Golay Filter')
#plt.show()
#plt.plot(x_loc, exp)
#plt.title('Exponential Smoothing')
#plt.show()

sav_exp = savgol_filter(exp, 100, 3)
exp_sav = exponential_smoothing(sav, 0.05)


c0 = [0.070248,
0.070248,
0.070248,
0.070248,
0.070248,
0.070248,
0.070248,
0.070248,
0.1405,
0.17562,
0.24587,
0.28099,
0.34836,
0.38053,
0.38053,
0.38055,
0.38056,
0.38059,
0.41726,
0.42182,
0.45521,
0.45541,
0.45573,
0.45622,
0.45677,
0.45759,
0.45898,
0.46037,
0.46273,
0.46736,
0.47483,
0.48755,
0.49809,
0.50679,
0.52229,
0.55966,
0.5925,
0.61698,
0.62152,
0.66246,
0.69983,
0.77458,
0.81195,
0.88669,
0.92406,
0.99881,
1.0362,
1.1109,
1.1483,
1.223,
1.2604,
1.3352,
1.3725,
1.4473,
1.4846,
1.5594,
1.5968,
1.6715,
1.7089,
1.7836,
1.821,
1.8957,
1.9331,
2.0078,
2.0452,
2.12,
2.1573,
2.2321,
2.2694,
2.3442,
2.3816,
2.4563,
2.4937,
2.5684,
2.6058,
2.6805,
2.7179,
2.7927,
2.83,
2.9048,
2.9421,
3.0169,
3.0543,
3.129,
3.1664,
3.2411,
3.2785,
3.3532,
3.3906,
3.4653,
3.5027,
3.5775,
3.6148,
3.6896,
3.7269,
3.8017,
3.8391,
3.9138,
3.9512,
4.0259,
4.0633,
4.138,
4.1754,
4.2501,
4.2875,
4.3623,
4.3996,
4.4744,
4.5118,
4.5714,
4.5957,
4.5994,
4.6016,
4.6033,
4.6051,
4.6063,
4.6392,
4.6442,
4.6855,
4.7226,
4.7982,
4.8372,
4.915,
4.954,
5.0319,
5.0708,
5.1487,
5.1876,
5.2655,
5.3045,
5.3824,
5.4213,
5.4992,
5.5382,
5.6161,
5.6551,
5.733,
5.772,
5.8499,
5.8889,
5.9668,
6.0058,
6.0838,
6.1227,
6.2007,
6.2397,
6.3177,
6.3566,
6.4346,
6.4736,
6.5516,
6.5906,
6.6686,
6.7076,
6.7856,
6.8246,
6.9026,
6.9416,
7.0196,
7.0586,
7.1366,
7.1756,
7.2536,
7.2927,
7.3707,
7.4097,
7.4877,
7.5268,
7.6048,
7.6438,
7.7219,
7.7609,
7.8389,
7.878,
7.956,
7.995,
8.0731,
8.1121,
8.1902,
8.2292,
8.3073,
8.3463,
8.4244,
8.4634,
8.5415,
8.5806,
8.6586,
8.6977,
8.7758,
8.8148,
8.8929,
8.9319,
9.01,
9.0491,
9.1272,
9.1662,
9.2443,
9.2834,
9.3615,
9.4005,
9.4786,
9.5177,
9.5958,
9.6348,
9.7129,
9.752,
9.8301,
9.8692,
9.9473,
9.9863,
10.064,
10.104,
10.182,
10.221,
10.299,
10.338,
10.416,
10.455,
10.533,
10.572,
10.65,
10.689,
10.768,
10.807,
10.885,
10.924,
11.002,
11.041,
11.119,
11.158,
11.236,
11.275,
11.353,
11.392,
11.471,
11.51,
11.588,
11.627,
11.705,
11.744,
11.822,
11.861,
11.939,
11.978,
12.057,
12.096,
12.174,
12.213,
12.291,
12.33,
12.408,
12.447,
12.525,
12.564,
12.642,
12.681,
12.76,
12.799,
12.877,
12.916,
12.994,
13.033,
13.111,
13.15,
13.228,
13.267,
13.345,
13.384,
13.462,
13.501,
13.58,
13.619,
13.697,
13.736,
13.814,
13.853,
13.931,
13.97,
14.048,
14.087,
14.165,
14.204,
14.282,
14.321,
14.399,
14.438,
14.516,
14.555,
14.633,
14.672,
14.75,
14.789,
14.867,
14.907,
14.985,
15.024,
15.102,
15.141,
15.219,
15.258,
15.336,
15.375,
15.453,
15.492,
15.57,
15.609,
15.687,
15.726,
15.804,
15.843,
15.921,
15.96,
16.038,
16.076,
16.154,
16.193,
16.271,
16.31,
16.388,
16.427,
16.505,
16.544,
16.622,
16.661,
16.739,
16.778,
16.856,
16.895,
16.973,
17.012,
17.09,
17.129,
17.206,
17.245,
17.323,
17.362,
17.44,
17.479,
17.557,
17.596,
17.674,
17.713,
17.79,
17.829,
17.907,
17.946,
18.024,
18.063,
18.141,
18.179,
18.257,
18.296,
18.373,
18.412,
18.49,
18.528,
18.606,
18.645,
18.722,
18.761,
18.838,
18.876,
18.954,
18.992,
19.069,
19.108,
19.185,
19.223,
19.3,
19.338,
19.415,
19.454,
19.53,
19.569,
19.645,
19.683,
19.76,
19.798,
19.874,
19.912,
19.989,
20.027,
20.103,
20.141,
20.216,
20.254,
20.33,
20.368,
20.444,
20.481,
20.557,
20.595,
20.67,
20.708,
20.783,
20.82,
20.895,
20.933,
21.008,
21.045,
21.12,
21.157,
21.232,
21.269,
21.343,
21.38,
21.454,
21.491,
21.565,
21.602,
21.676,
21.713,
21.787,
21.823,
21.895,
21.93,
21.999,
22.034,
22.103,
22.138,
22.207,
22.207,
22.207,
22.207,
22.207,
22.207,
22.207,
22.207]

M_loc = [1.397531452,
1.397532993,
1.397533574,
1.397534639,
1.397531055,
1.39752925,
1.397523354,
1.397524498,
1.397518059,
1.397524541,
1.397518059,
1.397537161,
1.397873944,
1.394108907,
1.393750817,
1.394660737,
1.393167819,
1.395688366,
1.429856034,
1.430304155,
1.466105786,
1.468948659,
1.470983663,
1.477590917,
1.46899763,
1.468635936,
1.466102656,
1.493449269,
1.443265719,
1.490111208,
1.483194916,
1.471401044,
1.465779866,
1.456082799,
1.453435187,
1.445056389,
1.448936898,
1.444422848,
1.453141043,
1.432948434,
1.438324784,
1.430957688,
1.430877668,
1.427376594,
1.428511974,
1.427033582,
1.426513749,
1.42887395,
1.428341402,
1.42939193,
1.429478149,
1.431485325,
1.430691173,
1.432441831,
1.432524869,
1.434187759,
1.433845818,
1.435859599,
1.43577006,
1.438052807,
1.437784553,
1.439721532,
1.439542614,
1.441914813,
1.441822242,
1.443666653,
1.443580046,
1.446310427,
1.446134206,
1.447444889,
1.44763028,
1.449744199,
1.449654441,
1.452214829,
1.452298704,
1.455208162,
1.455208114,
1.457160558,
1.457154642,
1.460075291,
1.460252531,
1.463347286,
1.463079658,
1.46689474,
1.467250127,
1.470001014,
1.470269398,
1.474082309,
1.474175739,
1.478177513,
1.478268097,
1.482187941,
1.482817903,
1.487454209,
1.488000029,
1.492737981,
1.492823628,
1.498290687,
1.499195904,
1.505477159,
1.506199149,
1.514201732,
1.514928894,
1.523139198,
1.524962797,
1.536010657,
1.538302505,
1.552577855,
1.558034259,
1.581557769,
1.529920018,
1.538171924,
1.5456315,
1.549515614,
1.556396659,
1.563155908,
1.684351368,
1.678539189,
1.866390874,
1.875718643,
1.834458613,
1.83568702,
1.71922905,
1.718717113,
1.681063393,
1.680756179,
1.660879115,
1.661390497,
1.649273623,
1.648762936,
1.64113816,
1.641547788,
1.635353498,
1.635250999,
1.630999586,
1.63151133,
1.628999259,
1.628180765,
1.626079761,
1.626693931,
1.624185232,
1.623878338,
1.62260075,
1.623010351,
1.6226575,
1.622042576,
1.621180764,
1.621796129,
1.621244467,
1.620732038,
1.620081098,
1.6204917,
1.620357322,
1.620151988,
1.620021001,
1.620431498,
1.620201813,
1.619790892,
1.620077974,
1.620077841,
1.620676949,
1.62067702,
1.620046403,
1.620251926,
1.620550006,
1.620447373,
1.621674828,
1.621571802,
1.621568846,
1.62167194,
1.621981501,
1.621981241,
1.622706145,
1.622500171,
1.622919922,
1.623126016,
1.623652471,
1.623652478,
1.624079711,
1.624182765,
1.624922834,
1.624716769,
1.625254538,
1.625357598,
1.626002108,
1.626002141,
1.626547154,
1.626650297,
1.627095933,
1.62699278,
1.627751508,
1.627648289,
1.627894832,
1.628101238,
1.629176627,
1.629073476,
1.629017589,
1.62912068,
1.630409811,
1.630100279,
1.630257771,
1.630567378,
1.632070195,
1.63186375,
1.631615307,
1.631718546,
1.632712299,
1.632815556,
1.633193423,
1.6328837,
1.63398772,
1.634710669,
1.635405288,
1.634269464,
1.634554246,
1.635690198,
1.63618503,
1.635255424,
1.634410839,
1.635237135,
1.639147201,
1.638527332,
1.636656514,
1.636553249,
1.641296557,
1.641399816,
1.636746555,
1.637056512,
1.640257073,
1.639740546,
1.639122408,
1.639432272,
1.63974709,
1.640057036,
1.642440985,
1.641821258,
1.642039278,
1.642039286,
1.641537549,
1.64205389,
1.641968569,
1.642175197,
1.643435727,
1.643229115,
1.644492754,
1.644492765,
1.6452432,
1.644830193,
1.644654487,
1.644860954,
1.64489496,
1.645204727,
1.645861422,
1.645861422,
1.646211517,
1.64610827,
1.646255148,
1.645945568,
1.647746887,
1.648056457,
1.649447657,
1.649447666,
1.648984684,
1.649191095,
1.649144087,
1.649247302,
1.650338179,
1.649822457,
1.649369269,
1.649472353,
1.650156971,
1.650053781,
1.650225666,
1.650328861,
1.650297684,
1.650400709,
1.651919047,
1.652434428,
1.652306259,
1.651893986,
1.652284235,
1.651975495,
1.652266322,
1.652369103,
1.652868777,
1.653177425,
1.65357724,
1.653165514,
1.65356824,
1.653877033,
1.653870972,
1.653665327,
1.653971354,
1.654074313,
1.653972054,
1.653971877,
1.654489765,
1.654592623,
1.654805409,
1.654702582,
1.654815373,
1.654918362,
1.655239839,
1.655137032,
1.655461569,
1.655564476,
1.655686635,
1.655584029,
1.655915022,
1.656017794,
1.656145961,
1.656249085,
1.656380812,
1.656380805,
1.656311192,
1.656310991,
1.656858823,
1.656756207,
1.657001021,
1.657000505,
1.65683744,
1.656735786,
1.656269583,
1.656474947,
1.657341784,
1.657238791,
1.656780572,
1.656780484,
1.657960299,
1.657755791,
1.656384189,
1.656486374,
1.656037846,
1.656037893,
1.657533528,
1.657533497,
1.656887775,
1.656887809,
1.657789062,
1.657890256,
1.660335668,
1.660030467,
1.6670453,
1.666741664,
1.67268969,
1.673702737,
1.67868164,
1.67827685,
1.683714274,
1.683915218,
1.688990129,
1.688789562,
1.694418423,
1.69441846,
1.698584627,
1.69798566,
1.701313906,
1.702010035,
1.706282405,
1.706381941,
1.71050395,
1.710209069,
1.713394773,
1.713791186,
1.715761541,
1.71566106,
1.716923066,
1.716628839,
1.719418662,
1.719808854,
1.721306797,
1.720823157,
1.724418702,
1.724420574,
1.724725053,
1.724442141,
1.725119486,
1.725312056,
1.725994443,
1.726365974,
1.726677616,
1.726959699,
1.725588245,
1.72513046,
1.721789367,
1.722061239,
1.721450932,
1.721359991,
1.718686848,
1.718591506,
1.717011318,
1.717460665,
1.713019687,
1.712134613,
1.704704161,
1.705589204,
1.700619459,
1.700278164,
1.692189959,
1.692452665,
1.685883021,
1.685632515,
1.673588488,
1.67383495,
1.65981207,
1.65981947,
1.644326706,
1.644188134,
1.625031207,
1.625184423,
1.600485275,
1.600485094,
1.576056348,
1.575776993,
1.53343699,
1.532984578,
1.459088378,
1.459031488,
1.387389772,
1.387288991,
1.321276862,
1.321290611,
1.321298504,
1.32125155,
1.321201407,
1.32111411,
1.321101746,
1.321066078]

theta = [0.3497575005,
0.3497575005,
0.3497575005,
0.3497575005,
0.3497575005,
0.3497575005,
0.3497575005,
0.3497575005,
0.3497575005,
0.3497575005,
0.3497575005,
0.3497575005,
0.3583611859,
0.3538480362,
0.3535389006,
0.3544130993,
0.3529527059,
0.3552769147,
0.2849922355,
0.284346237,
0.2181153382,
0.2179309738,
0.2179821854,
0.2177773425,
0.2180333976,
0.2180231551,
0.2178387944,
0.2180948529,
0.2178900049,
0.217971943,
0.217941216,
0.217941216,
0.2178900049,
0.2179514583,
0.217941216,
0.2179309738,
0.217941216,
0.2177671006,
0.2179002471,
0.217941216,
0.217941216,
0.2179309738,
0.217941216,
0.217941216,
0.2179309738,
0.2179309738,
0.217941216,
0.217941216,
0.2179309738,
0.217941216,
0.217941216,
0.2179309738,
0.217941216,
0.217941216,
0.2179309738,
0.2179309738,
0.217941216,
0.217941216,
0.2179309738,
0.217941216,
0.217941216,
0.217941216,
0.217941216,
0.217941216,
0.217941216,
0.217941216,
0.217941216,
0.217941216,
0.217941216,
0.2179309738,
0.217941216,
0.217941216,
0.217941216,
0.217941216,
0.217941216,
0.217941216,
0.217941216,
0.217941216,
0.217941216,
0.217941216,
0.217941216,
0.217941216,
0.217941216,
0.217941216,
0.217941216,
0.217941216,
0.217941216,
0.217941216,
0.217941216,
0.217941216,
0.217941216,
0.217941216,
0.217941216,
0.217941216,
0.217941216,
0.217941216,
0.217941216,
0.217941216,
0.217941216,
0.217941216,
0.217941216,
0.217941216,
0.217941216,
0.217941216,
0.217941216,
0.217941216,
0.217941216,
0.217941216,
0.217941216,
0.217941216,
0.2179309738,
0.217941216,
0.217941216,
0.2179514583,
0.2179309738,
0.2179309738,
0.1504064438,
0.1510436684,
0.08219050538,
0.07929507126,
0.08149317127,
0.08149016128,
0.07998224927,
0.07998325248,
0.07847451949,
0.07847451949,
0.07696395954,
0.07696395954,
0.07545558124,
0.07545457839,
0.07394536944,
0.07394536944,
0.07243733185,
0.07243632922,
0.07092845673,
0.07092945925,
0.06941773853,
0.06941773853,
0.06791118823,
0.06791018592,
0.06640078346,
0.06640178567,
0.06489153249,
0.06489053038,
0.06338242952,
0.06338242952,
0.06187347108,
0.06187347108,
0.06036365187,
0.0603646537,
0.05885497215,
0.05885497215,
0.05734642653,
0.05734542488,
0.05583600842,
0.05583700999,
0.05432772078,
0.05432772078,
0.05281855544,
0.05281955684,
0.0513095105,
0.05130850919,
0.0498005825,
0.0498005825,
0.04829176798,
0.04829176798,
0.04678306346,
0.04678306346,
0.04527146243,
0.04527346448,
0.04376396874,
0.04376396874,
0.04225557367,
0.04225457277,
0.04074527313,
0.04074527313,
0.03923606634,
0.03923706711,
0.03772694897,
0.03772794968,
0.03621791755,
0.0362189182,
0.03470896863,
0.03470896863,
0.03320009878,
0.03320109933,
0.03169230504,
0.03169130453,
0.030181582,
0.03018258245,
0.02867392909,
0.02867292868,
0.02716334028,
0.02716434064,
0.02565481411,
0.02565481411,
0.02414534605,
0.02414534605,
0.022635933,
0.02263693326,
0.02112757177,
0.02112757177,
0.01961825841,
0.01961825841,
0.01810998991,
0.01810998991,
0.01659976234,
0.01659976234,
0.01509157286,
0.01509257297,
0.01358141752,
0.01358041743,
0.01207329331,
0.01207329331,
0.01056319644,
0.01056319644,
0.009054623725,
0.009054923738,
0.007545571602,
0.007545071588,
0.006036536662,
0.006037036671,
0.004527415467,
0.004527115464,
0.003018204582,
0.003018104582,
0.001509100573,
0.001508900573,
0,
0,
0.001509200573,
0.001509000573,
0.003018204582,
0.003018104582,
0.004527015463,
0.00452671546,
0.006036536662,
0.006037036671,
0.007545571602,
0.00754517159,
0.009054623725,
0.009054923738,
0.01056319644,
0.0105641965,
0.01207229323,
0.01207229323,
0.01358241761,
0.01358141752,
0.01509057274,
0.01509157286,
0.0165987622,
0.01659976234,
0.01810998991,
0.01810998991,
0.01961825841,
0.01961825841,
0.02112657154,
0.02112757177,
0.02263693326,
0.02263793352,
0.02414534605,
0.02414434575,
0.02565381378,
0.02565381378,
0.02716334028,
0.02716434064,
0.02867492951,
0.02867392909,
0.030181582,
0.03018258245,
0.03169330554,
0.03169230504,
0.03320009878,
0.03320009878,
0.03470696743,
0.03470596683,
0.03622091952,
0.03622091952,
0.03772694897,
0.03772794968,
0.0392340648,
0.0392340648,
0.04074727479,
0.04074727479,
0.04225557367,
0.04225357188,
0.0437649697,
0.04376396874,
0.04527346448,
0.04527446551,
0.04678106127,
0.04678106127,
0.04829377031,
0.04829377031,
0.04979958126,
0.04979958126,
0.0513095105,
0.05130750787,
0.05281855544,
0.05281955684,
0.05432772078,
0.05432772078,
0.05583600842,
0.05583801155,
0.05734642653,
0.05734542488,
0.05885497215,
0.05885397041,
0.06036365187,
0.0603646537,
0.06187246917,
0.06187246917,
0.0633804255,
0.0633804255,
0.06489554093,
0.06489453882,
0.06639777684,
0.06639877905,
0.06791219054,
0.06791118823,
0.06941874094,
0.06941874094,
0.07092645169,
0.07092745421,
0.07243432396,
0.07243432396,
0.07394737492,
0.07394737492,
0.07545558124,
0.07545457839,
0.0769619536,
0.07696095063,
0.07847451949,
0.07847451949,
0.07998024286,
0.07998224927,
0.08149618126,
0.0814941746,
0.08398068105,
0.08398068105,
0.08737613772,
0.08737714155,
0.0946201253,
0.0946201253,
0.101876133,
0.101876133,
0.1091264614,
0.1091264614,
0.1163825538,
0.1163825538,
0.1236246532,
0.1236246532,
0.1308934482,
0.1308934482,
0.1381288209,
0.1381288209,
0.1453916918,
0.1453916918,
0.152631942,
0.152631942,
0.1598903944,
0.1598903944,
0.1671370706,
0.1671370706,
0.1743926177,
0.1743926177,
0.1816472842,
0.1816472842,
0.1889014495,
0.1889014495,
0.1961452984,
0.1961452984,
0.2033995891,
0.2033995891,
0.2106545241,
0.2106545241,
0.2178900049,
0.2178900049,
0.2251473545,
0.2251473545,
0.2324167873,
0.2324167873,
0.2396472786,
0.2396472786,
0.2469111927,
0.2469111927,
0.2541780989,
0.2541780989,
0.2613966404,
0.2613966404,
0.2686602649,
0.2686602649,
0.2759176583,
0.2759176583,
0.2831899957,
0.2831795809,
0.2904360081,
0.290425571,
0.2976663678,
0.2976663678,
0.3049337687,
0.3049442523,
0.3121757999,
0.312165292,
0.3194347973,
0.3194347973,
0.3266795906,
0.3266795906,
0.3339210424,
0.3339210424,
0.3411913436,
0.3412019554,
0.3484485441,
0.3484485441,
0.3586281584,
0.3586388379,
0.3498533024,
0.3498639472,
0.3498745921,
0.3498745921,
0.3498639472,
0.3498745921,
0.3498745921,
0.3498745921,
0.3498745921,
0.3498745921,
0.3498745921,
0.3498745921]

R = 4.287267 #ft
mu = np.arcsin(1/1.6)
#x_2 = [R*np.tan(90-np.radians(angle)) for angle in theta]

## Redefine nearfield_sig to be the exponential smoothed and savgol filtered data
nearfield_sig = exp_sav


## Calculates shock angle based on panel inclination and  local mach number
def theta_beta_mach(theta, M, gamma=1.4):
    # Define the equation to solve
    def equation(beta):
        left_hand_side = np.tan((theta))
        #right_hand_side = 2 * (1 / np.tan(beta)) * (
        #    (M**2 * (np.sin(beta)**2) - 1) /
        #    (M**2 * (gamma + np.cos(2 * beta)) + 2)
        #)
        right_hand_side = ((M**2 * np.sin(beta)**2 - 1)*(1/np.tan(beta))) / (((gamma+1)/2)*M**2 - M**2*np.sin(beta)**2 + 1)
        return left_hand_side - right_hand_side
    
    # Initial guess for beta (in radians)
    beta_initial_guess = (theta)  # A reasonable initial guess
    
    # Solve the equation
    beta_solution = fsolve(equation, beta_initial_guess)
    
    # Convert beta to degrees
    beta_degrees = beta_solution #np.degrees(beta_solution)
    
    return beta_degrees[0]

beta = []
for i in range(len(theta)):
    #if M_loc[i] < 1.6:
    #    beta.append(mu)
    #else:
    beta.append((theta_beta_mach(theta[i], M_loc[i])))
    
print("Beta: ", beta)
print("len_beta", len(beta))


## Defines new x start locatin from derived betas to match neirfield signature
x_2 = []
for i in range(len(c0)):
    x_2.append((R*np.tan((np.pi/2)-(beta[i])) + c0[i]))
#plt.plot(x_2,beta)
#plt.plot(x_2,theta)
#plt.show()

nodes_x = []
for value2 in x_2:
    closest_value = min(x_loc, key=lambda x:abs(x-value2))
    closest_index = x_loc.index(closest_value)
    nodes_x.append(closest_value)
    print(value2, closest_value, closest_index)

## Interpolate angles based on the nodes
angles_interpolated = (np.interp(x_loc, nodes_x , beta))
print("angles_interpolated: ", angles_interpolated)

## Calculate x_offset
x_offset = []
#R = 1.30675909 # meters
for i in range(len(angles_interpolated)):
    x_offset.append(abs((R/np.tan(angles_interpolated[i])) - (R/np.tan(mu))))
    if i >= len(angles_interpolated)/2:
        x_offset[i] = x_offset[i] / 3

## Calculate new x values
x_new = []
for i in range(len(x_offset)):
    x_new.append(x_loc[i] - x_offset[i])
    #print(x_loc[0][i], x_offset[i], x_new[i])


plt.plot(x_loc, nearfield_sig)
plt.plot(x_new, nearfield_sig)
plt.title("Original vs Corrected Nearfield Signature")
plt.legend(['Original', 'Corrected'])
plt.xlabel('X (ft)')
plt.ylabel('Pressure (Pa)')
#plt.xlim(0,75)
plt.show()
print(x_loc[-1] - x_loc[0], x_new[-1]-x_new[0])

p_new = []
for i in range(len(nearfield_sig)):
    p_new.append(nearfield_sig[i]*18823.1661 + 18823.1661)


# Read in a CSV file with 2 columns of data x,p
data = []
with open('studies/Goswift/results/Untitled spreadsheet - Sheet5.csv', 'r') as csvfile:
    reader = csv.reader(csvfile)
    next(reader)  # Skip the header row if it exists
    for row in reader:
        x = float(row[0])
        p = float(row[1])
        data.append((x, p))

plt.plot([x * 3.281 for x, _ in data], [p for _, p in data])
plt.plot(x_new, p_new)
plt.plot(x_loc, p_new)
plt.legend(['UNS3D', 'MachLine-Corrected', 'MachLine'])
plt.xlabel('X (ft)')
plt.ylabel('Pressure (Pa)')
plt.xlim(-5,80)
plt.show()