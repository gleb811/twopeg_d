#include <iomanip>
#include <string>
#include <stdio.h> 
#include <math.h>
#include <TLorentzVector.h>
#include <iostream>
#include "global.h"

using namespace std;


void fermi_bonn(Double_t R, Double_t R1, Double_t th_cos, Double_t ph) {



Double_t array[] = {
		0.,0.0000212,0.0001059,0.0002955,0.00063,0.0011487,
     0.0018873,0.0028801,0.0041582,0.0057497,0.0076795,0.0099690,
     0.0126359,0.0156946,0.0191559,0.0230270,0.0273117,0.0320108,
     0.0371220,0.0426400,0.0485571,0.0548634,0.0615465,0.0685926,
     0.0759860,0.0837102,0.0917470,0.1000778,0.1086834,0.1175439,
     0.1266394,0.1359500,0.1454558,0.1551369,0.1649743,0.1749491,
     0.1850426,0.1952373,0.2055162,0.2158629,0.2262616,0.2366973,
     0.2471563,0.2576250,0.2680910,0.2785424,0.2889683,0.2993584,
     0.3097032,0.3199937,0.3302220,0.3403805,0.3504618,0.3604606,
     0.3703704,0.3801863,0.3899037,0.3995183,0.4090266,0.4184252,
     0.4277112,0.4368822,0.4459359,0.4548708,0.4636848,0.4723770,
     0.4809465,0.4893924,0.4977144,0.5059122,0.5139856,0.5219348,
     0.5297601,0.5374620,0.5450408,0.5524971,0.5598319,0.5670460,
     0.5741401,0.5811155,0.5879735,0.5947151,0.6013415,0.6078541,
     0.6142542,0.6205431,0.6267222,0.6327930,0.6387572,0.6446159,
     0.6503709,0.6560237,0.6615756,0.6670282,0.6723834,0.6776425,
     0.6828066,0.6878780,0.6928575,0.6977472,0.7025485,0.7072628,
     0.7118914,0.7164361,0.7208984,0.7252796,0.7295811,0.7338039,
     0.7379501,0.7420208,0.7460173,0.7499411,0.7537935,0.7575758,
     0.7612892,0.7649351,0.7685146,0.7720289,0.7754796,0.7788674,
     0.7821937,0.7854598,0.7886667,0.7918154,0.7949073,0.7979431,
     0.8009241,0.8038512,0.8067257,0.8095480,0.8123197,0.8150414,
     0.8177143,0.8203391,0.8229172,0.8254489,0.8279356,0.8303778,
     0.8327764,0.8351323,0.8374463,0.8397193,0.8419515,0.8441445,
     0.8462986,0.8484146,0.8504933,0.8525352,0.8545414,0.8565121,
     0.8584484,0.8603507,0.8622199,0.8640565,0.8658607,0.8676338,
     0.8693760,0.8710880,0.8727704,0.8744236,0.8760484,0.8776451,
     0.8792143,0.8807565,0.8822724,0.8837624,0.8852267,0.8866661,
     0.8880810,0.8894718,0.8908390,0.8921831,0.8935046,0.8948038,
     0.8960810,0.8973367,0.8985714,0.8997853,0.9009790,0.9021528,
     0.9033068,0.9044417,0.9055577,0.9066554,0.9077349,0.9087966,
     0.9098409,0.9108681,0.9118784,0.9128722,0.9138499,0.9148115,
     0.9157574,0.9166883,0.9176040,0.9185050,0.9193915,0.9202638,
     0.9211218,0.9219661,0.9227971,0.9236149,0.9244196,0.9252118,
     0.9259914,0.9267586,0.9275138,0.9282570,0.9289889,0.9297093,
     0.9304183,0.9311165,0.9318039,0.9324808,0.9331472,0.9338033,
     0.9344492,0.9350855,0.9357119,0.9363288,0.9369365,0.9375350,
     0.9381244,0.9387048,0.9392768,0.9398401,0.9403951,0.9409419,
     0.9414804,0.9420111,0.9425339,0.9430493,0.9435568,0.9440572,
     0.9445504,0.9450365,0.9455154,0.9459873,0.9464527,0.9469115,
     0.9473636,0.9478095,0.9482490,0.9486823,0.9491096,0.9495308,
     0.9499463,0.9503559,0.9507599,0.9511585,0.9515517,0.9519394,
     0.9523219,0.9526993,0.9530715,0.9534388,0.9538010,0.9541582,
     0.9545110,0.9548590,0.9552026,0.9555416,0.9558762,0.9562063,
     0.9565322,0.9568540,0.9571718,0.9574854,0.9577950,0.9581008,
     0.9584028,0.9587010,0.9589953,0.9592862,0.9595734,0.9598572,
     0.9601374,0.9604142,0.9606878,0.9609581,0.9612251,0.9614889,
     0.9617496,0.9620072,0.9622619,0.9625136,0.9627624,0.9630083,
     0.9632515,0.9634919,0.9637294,0.9639644,0.9641968,0.9644266,
     0.9646538,0.9648787,0.9651010,0.9653208,0.9655384,0.9657536,
     0.9659665,0.9661772,0.9663856,0.9665918,0.9667960,0.9669980,
     0.9671980,0.9673960,0.9675919,0.9677858,0.9679778,0.9681678,
     0.9683560,0.9685424,0.9687270,0.9689097,0.9690907,0.9692701,
     0.9694476,0.9696236,0.9697979,0.9699705,0.9701415,0.9703110,
     0.9704790,0.9706454,0.9708104,0.9709740,0.9711361,0.9712968,
     0.9714561,0.9716139,0.9717705,0.9719256,0.9720795,0.9722321,
     0.9723834,0.9725335,0.9726824,0.9728299,0.9729764,0.9731216,
     0.9732658,0.9734088,0.9735506,0.9736913,0.9738310,0.9739696,
     0.9741072,0.9742437,0.9743792,0.9745136,0.9746472,0.9747797,
     0.9749112,0.9750418,0.9751714,0.9753002,0.9754280,0.9755550,
     0.9756811,0.9758063,0.9759307,0.9760542,0.9761769,0.9762987,
     0.9764197,0.9765400,0.9766594,0.9767781,0.9768960,0.9770132,
     0.9771297,0.9772454,0.9773604,0.9774746,0.9775882,0.9777011,
     0.9778132,0.9779248,0.9780356,0.9781458,0.9782553,0.9783642,
     0.9784725,0.9785801,0.9786871,0.9787936,0.9788994,0.9790046,
     0.9791092,0.9792133,0.9793168,0.9794197,0.9795221,0.9796239,
     0.9797252,0.9798259,0.9799261,0.9800258,0.9801250,0.9802237,
     0.9803218,0.9804194,0.9805166,0.9806132,0.9807094,0.9808051,
     0.9809003,0.9809951,0.9810894,0.9811832,0.9812766,0.9813695,
     0.9814619,0.9815540,0.9816456,0.9817368,0.9818275,0.9819178,
     0.9820077,0.9820971,0.9821862,0.9822749,0.9823632,0.9824510,
     0.9825385,0.9826255,0.9827122,0.9827985,0.9828844,0.9829699,
     0.9830551,0.9831399,0.9832243,0.9833083,0.9833920,0.9834753,
     0.9835583,0.9836409,0.9837232,0.9838051,0.9838866,0.9839678,
     0.9840487,0.9841293,0.9842095,0.9842893,0.9843689,0.9844481,
     0.9845269,0.9846055,0.9846838,0.9847617,0.9848393,0.9849166,
     0.9849936,0.9850702,0.9851466,0.9852226,0.9852984,0.9853738,
     0.9854489,0.9855238,0.9855983,0.9856726,0.9857466,0.9858202,
     0.9858935,0.9859666,0.9860394,0.9861119,0.9861841,0.9862560,
     0.9863276,0.9863990,0.9864701,0.9865409,0.9866114,0.9866816,
     0.9867516,0.9868213,0.9868907,0.9869599,0.9870288,0.9870974,
     0.9871658,0.9872339,0.9873017,0.9873693,0.9874366,0.9875036,
     0.9875704,0.9876369,0.9877032,0.9877692,0.9878349,0.9879004,
     0.9879656,0.9880306,0.9880953,0.9881598,0.9882240,0.9882880,
     0.9883517,0.9884152,0.9884784,0.9885414,0.9886041,0.9886666,
     0.9887289,0.9887909,0.9888526,0.9889141,0.9889754,0.9890364,
     0.9890972,0.9891578,0.9892181,0.9892782,0.9893381,0.9893977,
     0.9894571,0.9895163,0.9895751,0.9896338,0.9896923,0.9897505,
     0.9898084,0.9898663,0.9899238,0.9899811,0.9900381,0.9900949,
     0.9901515,0.9902079,0.9902640,0.9903199,0.9903756,0.9904311,
     0.9904863,0.9905413,0.9905961,0.9906507,0.9907051,0.9907592,
     0.9908131,0.9908667,0.9909202,0.9909735,0.9910265,0.9910793,
     0.9911320,0.9911844,0.9912365,0.9912885,0.9913402,0.9913917,
     0.9914429,0.9914941,0.9915450,0.9915957,0.9916462,0.9916965,
     0.9917465,0.9917963,0.9918459,0.9918953,0.9919445,0.9919935,
     0.9920423,0.9920908,0.9921392,0.9921873,0.9922352,0.9922829,
     0.9923304,0.9923778,0.9924250,0.9924718,0.9925184,0.9925650,
     0.9926113,0.9926574,0.9927033,0.9927490,0.9927945,0.9928398,
     0.9928848,0.9929298,0.9929744,0.9930189,0.9930631,0.9931072,
     0.9931511,0.9931948,0.9932384,0.9932817,0.9933248,0.9933676,
     0.9934103,0.9934528,0.9934951,0.9935372,0.9935791,0.9936209,
     0.9936625,0.9937039,0.9937450,0.9937860,0.9938269,0.9938675,
     0.9939079,0.9939480,0.9939880,0.9940279,0.9940675,0.9941070,
     0.9941463,0.9941854,0.9942243,0.9942630,0.9943016,0.9943399,
     0.9943781,0.9944161,0.9944538,0.9944914,0.9945288,0.9945660,
     0.9946032,0.9946401,0.9946767,0.9947133,0.9947497,0.9947858,
     0.9948218,0.9948576,0.9948932,0.9949287,0.9949640,0.9949992,
     0.9950341,0.9950688,0.9951034,0.9951378,0.9951720,0.9952059,
     0.9952399,0.9952736,0.9953071,0.9953405,0.9953737,0.9954067,
     0.9954395,0.9954723,0.9955047,0.9955371,0.9955693,0.9956013,
     0.9956331,0.9956648,0.9956963,0.9957277,0.9957588,0.9957898,
     0.9958206,0.9958513,0.9958817,0.9959121,0.9959424,0.9959725,
     0.9960024,0.9960321,0.9960617,0.9960911,0.9961203,0.9961494,
     0.9961783,0.9962071,0.9962356,0.9962641,0.9962924,0.9963205,
     0.9963484,0.9963762,0.9964039,0.9964314,0.9964587,0.9964859,
     0.9965129,0.9965398,0.9965665,0.9965932,0.9966196,0.9966459,
     0.9966722,0.9966981,0.9967239,0.9967497,0.9967753,0.9968007,
     0.9968260,0.9968511,0.9968760,0.9969009,0.9969255,0.9969500,
     0.9969745,0.9969987,0.9970227,0.9970466,0.9970703,0.9970939,
     0.9971175,0.9971409,0.9971641,0.9971873,0.9972102,0.9972330,
     0.9972556,0.9972782,0.9973007,0.9973229,0.9973450,0.9973670,
     0.9973890,0.9974107,0.9974323,0.9974539,0.9974753,0.9974965,
     0.9975176,0.9975386,0.9975594,0.9975801,0.9976006,0.9976211,
     0.9976414,0.9976615,0.9976816,0.9977015,0.9977214,0.9977410,
     0.9977606,0.9977801,0.9977994,0.9978185,0.9978376,0.9978566,
     0.9978753,0.9978941,0.9979127,0.9979311,0.9979495,0.9979677,
     0.9979858,0.9980038,0.9980217,0.9980395,0.9980571,0.9980747,
     0.9980921,0.9981093,0.9981265,0.9981436,0.9981605,0.9981773,
     0.9981941,0.9982107,0.9982272,0.9982436,0.9982600,0.9982761,
     0.9982921,0.9983081,0.9983239,0.9983395,0.9983550,0.9983705,
     0.9983858,0.9984011,0.9984163,0.9984314,0.9984464,0.9984612,
     0.9984760,0.9984906,0.9985052,0.9985196,0.9985340,0.9985482,
     0.9985623,0.9985763,0.9985904,0.9986042,0.9986179,0.9986316,
     0.9986451,0.9986586,0.9986719,0.9986852,0.9986984,0.9987115,
     0.9987245,0.9987373,0.9987501,0.9987627,0.9987753,0.9987879,
     0.9988003,0.9988126,0.9988248,0.9988369,0.9988490,0.9988610,
     0.9988729,0.9988846,0.9988963,0.9989079,0.9989194,0.9989308,
     0.9989421,0.9989533,0.9989645,0.9989756,0.9989866,0.9989974,
     0.9990081,0.9990188,0.9990295,0.9990401,0.9990506,0.9990610,
     0.9990713,0.9990815,0.9990916,0.9991017,0.9991117,0.9991217,
     0.9991316,0.9991413,0.9991511,0.9991607,0.9991702,0.9991797,
     0.9991891,0.9991985,0.9992077,0.9992168,0.9992259,0.9992349,
     0.9992439,0.9992527,0.9992614,0.9992702,0.9992789,0.9992875,
     0.9992961,0.9993045,0.9993128,0.9993211,0.9993293,0.9993375,
     0.9993456,0.9993536,0.9993616,0.9993694,0.9993773,0.9993851,
     0.9993927,0.9994004,0.9994080,0.9994155,0.9994230,0.9994304,
     0.9994377,0.9994450,0.9994522,0.9994593,0.9994664,0.9994733,
     0.9994803,0.9994872,0.9994940,0.9995008,0.9995075,0.9995142,
     0.9995207,0.9995273,0.9995337,0.9995402,0.9995466,0.9995529,
     0.9995592,0.9995653,0.9995714,0.9995775,0.9995835,0.9995894,
     0.9995953,0.9996012,0.9996070,0.9996127,0.9996184,0.9996241,
     0.9996296,0.9996352,0.9996406,0.9996461,0.9996514,0.9996568,
     0.9996621,0.9996673,0.9996725,0.9996776,0.9996827,0.9996877,
     0.9996926,0.9996975,0.9997025,0.9997073,0.9997121,0.9997169,
     0.9997216,0.9997262,0.9997308,0.9997353,0.9997398,0.9997443,
     0.9997487,0.9997531,0.9997575,0.9997618,0.9997661,0.9997703,
     0.9997745,0.9997787,0.9997827,0.9997868,0.9997908,0.9997948,
     0.9997987,0.9998026,0.9998065,0.9998103,0.9998141,0.9998178,
     0.9998214,0.9998251,0.9998287,0.9998323,0.9998358,0.9998393,
     0.9998428,0.9998462,0.9998496,0.9998530,0.9998564,0.9998596,
     0.9998628,0.9998661,0.9998693,0.9998724,0.9998756,0.9998786,
     0.9998817,0.9998847,0.9998877,0.9998907,0.9998937,0.9998965,
     0.9998993,0.9999022,0.9999050,0.9999078,0.9999105,0.9999133,
     0.9999159,0.9999186,0.9999212,0.9999238,0.9999263,0.9999288,
     0.9999313,0.9999338,0.9999363,0.9999387,0.9999411,0.9999435,
     0.9999459,0.9999482,0.9999505,0.9999527,0.9999549,0.9999571,
     0.9999593,0.9999615,0.9999636,0.9999658,0.9999679	};


Double_t theta,phi, momentum;

for (Int_t j=0; j<=1000; j++) {
if ((R >= array[j])&&(R < array[j+1])) momentum = 0.001*double(j)+0.001*R1;
};

//theta =  acos(1.-2.*((double) rand()/(RAND_MAX)));
theta = acos(th_cos);

phi = ph;

px_fermi = momentum*sin(theta)*cos(phi);
py_fermi = momentum*sin(theta)*sin(phi);
pz_fermi = momentum*cos(theta);

};
