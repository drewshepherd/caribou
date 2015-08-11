/*FissionHeatMaterial source file
	*
	*inputs average burnup, enrichment and pellet radius and interpolates to give the three flux parameter constants (kappa, beta and lambda) in the fission 		heat equation (Q_fission), then it outputs Q_fission
	*
	*written by Drew Shepherd and Kyle Gamble
*/

#include "FissionHeatMaterial.h"
#include <math.h>
#include "Function.h"

template<>
InputParameters validParams<FissionHeatMaterial>()
{
  InputParameters params = validParams<Material>();
	params.addRequiredParam<FunctionName>("linear_power", "The linear element power function (W/m)");
	params.addRequiredParam<bool>("model_Qfission", "Set true to calculate the heat produced by fission");
	params.addRequiredCoupledVar("burnup_avg","the average burnup over the radial surface of the pellet at an instant in time");
	params.addRequiredParam<Real>("enrichment", "The percentage enrichment of the fuel");	
	params.addRequiredParam<Real>("pellet_radius","Radius of the pellet");
	params.addRequiredParam<Real>("ratio","The ratio of thermal power to fission power, default = 0.925");
	params.addRequiredParam<Real>("initial_fuel_density", "The initial density of the fuel(kg/m^3)");
	params.addRequiredParam<Real>("initial_qfission", "Initial heating term");
	params.addRequiredParam<Real>("initial_fuel_area", "Initial area of the fuel");
	params.addRequiredParam<bool>("is_3D", "Is the geometry in 2 or 3 dimensions, true if it is 3D, false if it is 2D");
	params.addRequiredParam<bool>("model_plate_fuel", "Is the geometry in plate fuel or a rod");
	return params;
}
//Constructor - obtains the values of burnup, enrichment and pellet radius from the input file
FissionHeatMaterial::FissionHeatMaterial(const InputParameters & parameters)
  :Material(parameters),

  _linear_power(&getFunction("linear_power")),
	_model_Qfission(getParam<bool>("model_Qfission")),
	_burnup_avg(coupledValue("burnup_avg")),
	_enrichment_property(getParam<Real>("enrichment")),
	_pellet_radius_property(getParam<Real>("pellet_radius")),  
	_ratioProperty(getParam<Real>("ratio")),  
	_initial_density(getParam<Real>("initial_fuel_density")),
	_initial_qfission(getParam<Real>("initial_qfission")),
	_area_fuel(getParam<Real>("initial_fuel_area")),
	_is_3D(getParam<bool>("is_3D")),
	_model_plate_fuel(getParam<bool>("model_plate_fuel")),

	_density(getMaterialProperty<Real>("density")),
	_q_fission(declareProperty<Real>("q_fission")),
	_q_fission_old(declarePropertyOld<Real>("q_fission")),
	_ratio(declareProperty<Real>("_ratio")),
	_pellet_radius(declareProperty<Real>("pellet_rad")),
	_enrichment(declareProperty<Real>("enrich"))
{
}

//Interpolates the burnup and enrichment values, this function is called in Q_fission() and used to find the appropriate beta, kappa and lambda constants
double FissionHeatMaterial::interp2(const double cons_array[][10])
{
	//defining arrays and their size
	const int enrichment_size = 10;
	const int burnup_size = 21;	
	const double enrichment_array[enrichment_size] = {0.71, 1, 1.5, 2, 2.5, 3, 3.5, 4, 5, 6};
	const int burnup_array[burnup_size] = {0,	48,	96,	144, 192,	240,	288,	336,	384,	432,	480,	528,	576,	624,	672,	720,	768,	816,	864,	912,	960};

	int i = 1;
	double enrich1;
	double enrich2;
	int burnup1;
	int burnup2;

	int enrich_index1;
	int enrich_index2;
	int bu_index1;
	int bu_index2;

	//determining the enrichment coordinates
	for (i = 1; i < enrichment_size; ++i)
	{
		if (_enrichment_property < enrichment_array[i])
		{
			//Enrichment values
			enrich_index1 = i-1;
			enrich_index2 = i;
			enrich1 = enrichment_array[enrich_index1];
			enrich2 = enrichment_array[enrich_index2]; 
			break;
		}
		else if (_enrichment_property == enrichment_array[9])
		{
			enrich_index1 = i;
			enrich_index2 = i;
			enrich1 = enrichment_array[9];
			enrich2 = enrichment_array[9];
			break;
		}
	}

	const Real burnup_avg_ = _burnup_avg[_qp];

	//determining the burnup coordinates
	for (i = 1; i < burnup_size; ++i)
	{
		if (burnup_avg_ < burnup_array[i])
		{	
			//Burnup values
			bu_index1 = i-1;
			bu_index2 = i;
			burnup1 = burnup_array[bu_index1];
			burnup2 = burnup_array[bu_index2];
			break;		
		}
		else if (burnup_avg_ == burnup_array[20])
		{
			bu_index1 = i;
			bu_index2 = i;
			burnup1 = burnup_array[20];
			burnup2 = burnup_array[20];
			break;		
		}
	}

	//determining the constant's value using enrichment as a variable - this is the first interpolation
	const double cons_lowburnup = (cons_array[bu_index1][enrich_index2] - cons_array[bu_index1][enrich_index1]) / (enrich2 - enrich1) * (_enrichment_property - enrich1) + cons_array[bu_index1][enrich_index1];
	//the first interpolation done again
	const double cons_highburnup = (cons_array[bu_index2][enrich_index2] - cons_array[bu_index2][enrich_index1]) / (enrich2 - enrich1) * (_enrichment_property - enrich1) + cons_array[bu_index2][enrich_index1];
	
	//determining the constant's value using burnup as a variable - this is the second interpolation
	const double _greek_cons = (cons_highburnup - cons_lowburnup) / (burnup2 - burnup1) * (burnup_avg_ - burnup1) + cons_lowburnup;

return _greek_cons; //passes this value to findConstants();
}

//Taken from HORSE - QFission, this isone of the bessel functions
Real
FissionHeatMaterial::bessel0(double x)
{
	double ax, ans;
	double y;  //Accumulate polynomials in double precision.
	
	if ((ax = fabs(x)) < 3.75)  // Polynomial fit
	{
		y = x / 3.75;
		y*=y;
		ans = 1.0 + y * (3.5156229 + y * (3.0899424 + y * (1.2067492 + y * (0.2659732 + y * (0.360768e-1 + y * 0.45813e-2)))));
	}
	else
	{
		y = 3.75/ax;
		ans = (std::exp(ax) / std::sqrt(ax)) * (0.39894228 + y * (0.1328592e-1 + y * (0.225319e-2 + y * (-0.157565e-2 + y * (0.916281e-2 + y * (-0.2057706e-1 + y * (0.2635537e-1 + y * (-0.1647633e-1 + y * (0.392377e-2)))))))));
	}
	return ans;
}

//Taken from HORSE - QFission, this is one of the bessel functions
Real 
FissionHeatMaterial::bessel1(double x)
{
	double ax, ans;
	double y;	//Accumulate polynomials in double precision.
	
	if ((ax = fabs(x)) < 3.75) 	//Polynomial fit
	{
		y = x / 3.75;
		y*=y;
		ans = ax*(0.5 + y * (0.87890594 + y * (0.51498869 + y * (0.15084934 + y * (0.2658733e-1 + y * (0.301532e-2 + y * 0.32411e-3))))));
	}
	else
	{
		y = 3.75 / ax;
		y*=y;
		ans = 0.2282967e-1 + y * (-0.2895312e-1 + y * (0.1787654e-1 - y * 0.420059e-2));
		ans = 0.39894228 + y * (-0.3988024e-1 + y * (-0.362018e-2 + y * (0.163801e-2 + y * (-0.1031555e-1 + y * ans))));
		ans *= (std::exp(ax) / std::sqrt(ax));
	}
	return x < 0.0 ? -ans : ans;
}

//calculates the values of kappa, beta and lambda
Real FissionHeatMaterial::Q_fission()
{
//these are the matrices of the flux parameters, enrichment, burnup and pellet_radius are used to interpolate the correct values. The ten enrichment values are the columns, the 21 burnup values are the rows. Each of these nine matrices corresponds to a different pellet_radius
	const double r4_beta[21][10] = {{0,0,0,0,0,0,0,0,0,0},
{0.1697,0.1324,0.0947,0.0755,0.658,0.615,0.0614,0.0629,0.07,0.0817},
{0.287,0.2383,0.1775,0.1402,0.116,0.1024,0.0921,0.0874,0.0864,0.0912},
{0.371,0.3236,0.2502,0.1987,0.166,0.1428,0.1276,0.1159,0.1054,0.1044},
{0.4293,0.3933,0.3155,0.2596,0.2154,0.1842,0.1628,0.145,0.1267,0.1199},
{0.4695,0.4471,0.3766,0.3098,0.2608,0.2239,0.1969,0.1747,0.1495,0.1361},
{0.497,0.4857,0.4296,0.3613,0.3039,0.2611,0.2317,0.2045,0.1725,0.1539},
{0.5331,0.5124,0.4743,0.4089,0.3474,0.2993,0.264,0.2323,0.1956,0.1726},
{0.5233,0.259,0.5096,0.4516,0.3893,0.3372,0.2969,0.2627,0.2177,0.1908},
{0.5278,0.537,0.5358,0.4921,0.4296,0.3733,0.3279,0.2904,0.2412,0.2089},
{0.5319,0.541,0.5519,0.526,0.4673,0.4093,0.3599,0.3175,0.2637,0.2274},
{0.5333,0.5411,0.5578,0.5515,0.5031,0.4442,0.3925,0.345,0.2852,0.2459},
{0.5333,0.5406,0.5606,0.5671,0.5349,0.4784,0.4237,0.3737,0.308,0.2642},
{0.5346,0.5403,0.5584,0.5766,0.5604,0.5117,0.4553,0.4009,0.3306,0.2839},
{0.5344,0.5385,0.5553,0.5779,0.5804,0.5393,0.4862,0.43,0.3522,0.3016},
{0.5333,0.5382,0.5515,0.575,0.5901,0.5658,0.5159,0.458,0.3741,0.3186},
{0.5325,0.5359,0.5469,0.5692,0.5934,0.5867,0.544,0.4868,0.3967,0.3377},
{0.5332,0.5334,0.542,0.5616,0.5889,0.5986,0.5697,0.5143,0.4201,0.356},
{0.5308,0.532,0.5397,0.555,0.582,0.6027,0.5897,0.5411,0.4645,0.3744},
{0.5289,0.5308,0.5357,0.5482,0.5722,0.6016,0.6041,0.5653,0.5069,0.3916},
{0.5285,0.5287,0.5319,0.5421,0.5628,0.5939,0.6113,0.5868,0.5463,0.4392}};
	const double r4_kappa[21][10] = {{107.75,120.86,140.65,157.97,174.12,188.87,202.86,215.97,240.66,263.15},
{113.57,118.49,131.23,145.19,158.54,171.18,183.4,194.67,215.69,234.67},
{118.95,119.49,127.66,139.67,152.62,165.42,177.66,189.21,210.97,230.59},
{122.91,119.76,123.36,133.14,145.56,158.05,170.5,182.6,204.86,225.34},
{127.1,121.24,119.55,126.74,137.75,150.1,162.66,174.78,197.99,219.2},
{131.99,123.93,116.96,120.48,129.64,141.46,153.96,166.34,190.34,212.36},
{136.42,128.1,116.18,115.14,121.76,132.66,144.86,157.38,182.09,204.95},
{140.2,132.81,117.36,110.7,114.05,123.45,135.11,147.78,173.27,197.09},
{143.2,137.13,120.66,108.38,107.19,114.35,125.27,137.86,163.94,188.59},
{145.12,140.83,125.58,108.9,101.71,105.65,114.96,127.22,154.02,176.58},
{146.55,143.63,130.79,111.63,98.02,97.32,104.48,115.95,143.44,170.06},
{147.21,145.38,135.7,116.55,97.4,90.25,94.15,104.42,132.25,159.89},
{147.59,146.65,139.89,122.59,99.56,85.31,84.35,92.36,120.3,149.25},
{147.75,147.29,143.04,128.95,104.86,83.5,75.4,79.76,107.47,137.75},
{147.78,147.61,145.22,134.77,111.89,85.37,68.56,67.14,93.44,125.26},
{147.59,147.74,146.42,139.29,119.75,90.89,65.49,54.78,77.65,111.75},
{147.35,147.6,147.08,142.6,127.65,98.96,66.7,44,59.27,96.93},
{147.13,147.41,147.43,144.77,134.06,108.68,72.9,36.71,35.16,79.76},
{146.88,147.2,147.56,146.21,138.97,118.43,82.89,37.51,35.16,58.71},
{146.48,146.87,147.36,146.87,142.35,127.15,94.68,46.59,35.16,26.07},
{146.3,146.54,147.13,147.16,144.56,134.09,107.06,60.72,35.16,26.07}};
	const double r4_lambda[21][10] = {{100,100,100,100,100,100,100,100,100,100},
{9726.5,9253.35,8138.53,6959.92,5846.6,4964.49,4428.31,3960.19,3374.88,3039.43},
{10168.3,10009.69,9519.63,8760.13,7907.99,7071.8,6206.79,5523.93,4571.25,3904.66},
{10312.96,10315.87,10021.77,9534.8,9085.3,8290.38,7598.38,6892.21,5660.8,4810.13},
{10314.25,10427.46,10239.94,10023.85,9672.5,9138.86,8587.62,7896.15,6693.76,5703.59},
{10326.87,10448.67,10416.46,10294.23,10064.78,9706.33,9256.58,8687.08,7584.52,6533.15},
{10329.9,10423.77,10497.4,10494.56,10349.84,10091.62,9762,9331.82,8358,7311.94},
{10311.11,10405.2,10514.95,10597.54,10531.57,10367.38,10109.45,9770.02,8997.33,8081.02},
{10302.01,10382.68,10491.84,10624.16,10664.7,10576.81,10421.96,10196.63,9521.23,8718.21},
{10262.09,10343.42,10473.49,10653.89,10731.52,10727.16,10602.51,10459.92,9999.68,9252.3},
{10267.15,10318.41,10433.59,10631.11,10751.32,10803.55,10756.68,10652.75,10361.18,9746.05},
{10251.53,10274.91,10361.55,10584.95,10738.28,10845.24,10869.53,10812.38,10632.71,10176.65},
{10243.76,10264.48,10334.91,10502.59,10686.92,10835.07,10918.21,10936.26,10878.85,10567.79},
{10244.89,10257.26,10307.38,10452.33,10636.2,10834.08,10941.19,10995.96,11056.11,10893.65},
{10262.43,10256.36,10295.47,10395.03,10582.13,10759.84,10924.3,11043.88,11195.75,11109.29},
{10243.66,10265.88,10276.36,10352.97,10503.46,10708.81,10902.14,11053.17,11267.03,11286.33},
{10246.55,10250.25,10257.46,10321.09,10457.89,10636.62,10831.67,11039.09,11319.11,11452.5},
{10266.74,10246.16,10244.26,10283.08,10389.57,10560.68,10764.02,10975.86,11346.6,11569.59},
{10255.39,10245.06,10253.41,10278.99,10343.39,10487.63,10682.49,10930.92,12689.21,11657.66},
{10241.48,10254.98,10244.31,10253.46,10292.81,10428.7,10597.3,10846.98,13767.27,11670.47},
{10263.9,10243.58,10234.37,10238.18,10262.9,10378.15,10530.33,10751.49,14491.31,13741.15}};
	const double r6_beta[21][10] = {{0,0,0,0,0,0,0,0,0,0},
{0.3731,0.2884,0.2063,0.1611,0.1369,0.1243,0.119,0.1182,0.1287,0.1479},
{0.6247,0.5169,0.3865,0.3068,0.255,0.2215,0.1992,0.1837,0.1719,0.1775},
{0.8058,0.6967,0.542,0.438,0.3662,0.3176,0.281,0.2543,0.2242,0.214},
{0.9394,0.8436,0.6782,0.5581,0.4706,0.4078,0.3633,0.324,0.2789,0.256},
{1.0341,0.9595,0.8035,0.6683,0.5665,0.4978,0.4394,0.391,0.3338,0.2996},
{1.1049,1.0525,0.9095,0.7687,0.6579,0.4931,0.5115,0.4558,0.3895,0.3454},
{1.1523,1.1226,1.0043,0.8836,0.7426,0.574,0.58,0.5194,0.4408,0.3898},
{1.1831,1.1721,1.0843,0.9492,0.8242,0.6504,0.6463,0.5788,0.4918,0.4332},
{1.2038,1.204,1.1514,1.0286,0.9012,0.7229,0.7094,0.6359,0.541,0.4755},
{1.2173,1.2254,1.2001,1.1014,0.9782,0.7928,0.7716,0.6919,0.5877,0.5185},
{1.2226,1.2366,1.234,1.1562,1.0439,0.8598,0.8295,0.7439,0.6317,0.558},
{1.2281,1.2404,1.2559,1.2124,1.1059,0.9274,0.8871,0.7962,0.677,0.595},
{1.2286,1.2404,1.2657,1.2504,1.1663,0.989,0.9437,0.8468,0.7187,0.6334},
{1.2291,1.2395,1.2673,1.2737,1.2145,1.05,0.9998,0.8961,0.7602,0.6684},
{1.228,1.2358,1.2649,1.2873,1.2538,1.1079,1.0534,0.9449,0.8003,0.7038},
{1.2243,1.2321,1.2585,1.2894,1.283,1.1606,1.1054,0.9933,0.8396,0.738},
{1.2209,1.2284,1.2518,1.2856,1.2992,1.21,1.1544,1.0419,0.8783,0.7713},
{1.2169,1.2237,1.2427,1.2784,1.3058,1.2492,1.2004,1.087,0.9162,0.8025},
{1.2141,1.2165,1.2358,1.2673,1.3023,1.282,1.2417,1.1328,0.9549,0.8342},
{1.2081,1.213,1.229,1.2544,1.2953,1.3041,1.2753,1.1748,0.9943,0.8654}};
	const double r6_kappa[21][10] = {{89.94,100.71,117,131.13,143.93,155.8,166.8,177.12,196.33,213.64},
{98.13,101.64,111.62,122.61,133.17,143.06,152.2,160.75,176.39,190.47},
{103.93,104.05,110.08,119.53,129.64,139.56,149.06,157.89,174.05,188.62},
{108.1,105.52,108.07,115.64,124.99,134.81,144.37,153.49,170.46,185.57},
{111.9,107.34,106.15,111.64,119.93,129.39,139.11,148.33,165.92,181.7},
{115.56,109.75,105.07,107.66,114.81,123.78,133.37,142.68,160.77,177.26},
{119.14,112.7,104.7,104.37,109.76,117.93,127.26,136.67,155.27,172.35},
{122.17,116.1,105.49,102.11,104.82,112.05,120.97,130.36,149.27,166.99},
{124.47,119.24,107.38,100.03,100.51,106.22,114.6,123.78,143.01,161.3},
{126.23,122.01,110.19,99.61,96.92,100.61,108.12,113.93,136.41,155.24},
{127.5,124.34,113.38,100.63,94.25,95.32,101.63,109.87,129.59,148.98},
{128.23,126.03,112.79,103.24,92.92,90.8,95.27,102.69,122.46,142.29},
{128.75,127.23,119.9,106.07,92.93,87.23,98.2,95.39,114.97,135.23},
{128.94,128.05,122.54,109.83,94.65,84.88,83.66,88.19,107.09,129.8},
{129.05,128.54,124.68,113.76,97.57,84.12,79.08,81.06,98.9,120.32},
{129.05,128.77,126.18,117.49,101.57,85,75.7,74.41,90.33,112.21},
{128.88,128.85,127.23,120.6,106.17,87.7,73.93,68.39,81.39,103.64},
{128.73,128.87,127.91,123.13,110.74,91.88,74.18,63.65,71.97,94.52},
{128.55,128.73,128.24,125.02,115.13,97.03,76.32,60.55,62.13,84.73},
{128.34,128.54,128.47,126.3,118.75,102.57,80.57,59.83,51.8,74.17},
{128.1,128.35,128.51,127.15,121.75,108.03,86.17,61.67,40.93,62.38}};
	const double r6_lambda[21][10] = {{100,100,100,100,100,100,100,100,100,100},
{8035.63,7627.75,6773.99,5774.12,4854.93,4099.27,3532.27,3079.03,2558.4,2290.82},
{8391.87,8249.98,7833.6,7279.66,6605.48,5910.32,5248.69,4629.28,3680.14,3114.5},
{8514.39,8454.69,8256.35,7932.64,7470.87,6954.79,6368,5793.5,4757.82,3954.23},
{8552.43,8549.62,8470.89,8289.02,7980.35,7597.25,7169.38,6657.44,5662.16,4769.43},
{8532.3,8570.36,8612.79,8495.28,8298.04,8033.49,7711.45,7301.92,6415.82,5535.1},
{8524.06,8568.31,8639.27,8607.28,8516.34,8327.27,8094.11,7783.13,7060.32,6229.98},
{8506.05,8557.83,8646.66,8929.32,8644.33,8540.97,8366.76,8170.88,7548.82,6824.12},
{8481.78,8532.4,8638.23,8694.44,8738.04,8683.24,8591.38,8452.5,7980.38,7350.59},
{8469.16,8508.65,8622.09,8697.32,8780.95,8768.9,8538.18,8660.61,8303.47,7794.45},
{8461.99,8494.09,8580.23,8685.51,8818.64,8821.82,8850.38,8819.66,8579.16,8191.27},
{8445.24,8474.16,8548.36,8666.59,8785.7,8854.74,8907.64,8920.78,8791.17,8506.07},
{8445.15,8459.49,8522.61,8622.83,8745.94,8851.32,8938.92,8997.99,8974.61,8754.19},
{8434.82,8445.32,8495.44,8587.06,8723.39,8832.36,8943.08,9045.49,9095.87,9004.51},
{8437.16,8439.31,8471.17,8548.94,8672.28,8803.69,8938.44,9057.25,9192.84,9181.62},
{8439.19,8432.65,8456.3,8520.1,8627.35,8756.17,8912.72,9051.68,9286.44,9330.62},
{8431.15,8427.56,8443.52,8485.61,8584.87,8709.25,8869.38,9028.8,9286.44,9451.65},
{8424.98,8425.32,8436.18,8463.18,8537.87,8657.58,8820.15,8989.95,9286.74,9526.09},
{8422.47,8422.12,8423.19,8450.73,8508.76,8613.19,8758.8,8933.06,9277.81,9559.71},
{8419.65,8409.88,8420.23,8433.19,8471.8,8565.5,8704,8874.21,9250.34,9588.84},
{8411.84,8412.42,8419.93,8417.77,8456.09,8519.68,8643.55,8807.27,9204.47,9594.29}};
	const double r9_beta[21][10] = {{0,0,0,0,0,0,0,0,0,0},
{0.7934,0.6239,0.4502,0.3526,0.2927,0.2595,0.2408,0.2332,0.2384,0.2605},
{1.2995,1.0882,0.8319,0.6715,0.5635,0.4894,0.4355,0.3983,0.3563,0.3464},
{1.6555,1.4414,1.1446,0.9433,0.8053,0.7047,0.628,0.5688,0.4908,0.4494},
{1.9188,1.7202,1.4117,1.1819,1.0182,0.8995,0.8052,0.7314,0.6255,0.5601},
{2.1223,1.945,1.6345,1.3909,1.2101,1.0745,0.9665,0.8821,0.7579,0.6694},
{2.2798,2.1281,1.8323,1.5767,1.3822,1.2324,1.162,1.0216,0.8811,0.7793},
{2.3983,2.275,2.0035,1.7459,1.5387,1.3814,1.2543,1.1534,0.999,0.885},
{2.4869,2.3932,2.1527,1.8957,1.6837,1.5152,1.3794,1.2721,1.1077,0.9838},
{2.5486,2.4849,2.2812,2.0338,1.8155,1.6398,1.4976,1.3839,1.2084,1.0796},
{2.5936,2.5533,2.3897,2.1577,1.937,1.755,1.6051,1.4875,1.3047,1.1705},
{2.627,2.6034,2.4775,2.2662,2.0501,1.8628,1.7083,1.5842,1.3942,1.2511},
{2.6476,2.6359,2.5466,2.364,1.1532,1.9613,1.8052,1.6747,1.4786,1.3306},
{2.6593,2.6581,2.6012,2.4471,2.2475,2.0545,1.8919,1.7573,1.5527,1.4066},
{2.6635,2.6688,2.6384,2.5177,2.331,2.1411,1.9747,1.837,1.6289,1.474},
{2.6619,2.6735,2.6641,2.5745,2.4092,2.2204,2.0544,1.9123,1.6938,1.5388},
{2.6592,2.6733,2.6768,2.6174,2.4748,2.2959,2.1247,1.9805,1.7581,1.6},
{2.6522,2.666,2.6854,2.6489,2.532,2.3652,2.1942,2.0444,1.8194,1.6556},
{2.6433,2.658,2.6801,2.6681,2.5777,2.4248,2.2607,2.109,1.8759,1.709},
{2.632,2.2465,2.6753,2.6796,2.6147,2.4799,2.3219,2.1675,1.9276,1.7572},
{2.6192,2.6318,2.6651,2.6794,2.6407,2.5282,2.3775,2.2232,1.9765,1.8039}};
	const double r9_kappa[21][10] = {{75.3,83.99,96.81,107.88,117.66,126.47,134.65,142.13,155.76,167.8},
{86.19,88.43,95.56,103.45,111.06,118.06,124.45,130.23,140.47,149.13},
{92.96,92.48,96.59,103.11,110.17,117.11,123.66,129.72,140.44,149.57},
{96.97,95.19,96.85,101.95,108.36,115.09,121.63,127.82,139.09,148.72},
{100.33,97.42,96.92,100.6,106.24,112.6,119.05,125.35,136.96,147.08},
{103.19,99.49,96.96,99.2,104,109.92,116.24,122.53,134.43,144.95},
{105.6,101.48,97.27,97.91,101.74,107.23,113.29,119.51,135.7,142.47},
{107.69,103.37,97.79,96.95,99.55,104.45,110.2,116.35,128.52,139.72},
{109.3,105.21,98.76,96.18,97.6,101.67,107.07,113.09,125.3,136.77},
{110.5,106.89,99.91,95.74,95.82,98.95,103.96,109.73,121.94,133.61},
{111.54,108.35,101.21,95.71,94.34,96.38,100.84,106.3,118.4,130.37},
{112.34,109.63,102.62,96.06,93.19,94.09,97.81,102.9,114.85,126.9},
{112.91,110.66,104.13,96.84,92.43,91.98,94.89,99.46,111.14,123.36},
{113.29,111.52,105.6,97.87,92.16,90.24,92.04,96.06,107.35,119.71},
{113.54,112.16,106.95,99.21,92.29,88.89,89.48,92.77,103.52,115.86},
{113.73,112.65,108.23,100.66,92.87,87.98,87.25,89.57,99.62,111.98},
{113.82,113.01,109.33,102.28,93.89,87.6,85.29,86.59,95.39,107.96},
{113.85,113.24,110.26,103.85,95.22,87.7,83.83,83.85,91.76,103.86},
{113.83,113.41,111,105.34,96.78,88.26,82.88,81.44,87.81,99.65},
{113.77,113.49,111.62,106.76,98.56,89.32,82.43,79.45,83.93,95.32},
{113.69,113.5,112.08,107.96,100.36,90.81,82.55,77.89,80.17,90.95}};
	const double r9_lambda[21][10] = {{100,100,100,100,100,100,100,100,100,100},
{6687.88,6384.86,5712.99,4930.77,4159.78,3525.63,3003.87,2608.6,2105.1,1814.33},
{6967.77,6852.88,6541.04,6107.31,5591.36,5047.73,4509.4,4004.95,3183.05,2622.58},
{7061.32,7022.76,6855.46,6594.78,6261.09,5047.73,5432.03,7974.75,4127.56,3413.47},
{7086.31,7086.62,7023.27,6864.45,6636.44,5832.6,6021.65,5648.91,4872.99,4133.02},
{7087.33,7113.41,7087.67,6864.45,6887.6,6362.34,6424.94,6137.59,5474.74,4753.97},
{7091.05,7114.44,7127.24,7014.52,7023.76,6677.96,6718.23,6488.1,5940.28,5296.25},
{7094.3,7104.64,7141.8,7093.98,7117.96,6898.68,6929.16,6763.98,6318.72,5757.13},
{7076.15,7095.11,7141.8,7149.71,7176.85,7056.48,7070.84,6967.58,6628.01,6136.47},
{7056.57,7081.28,7136.84,7173.43,7212.78,7158.74,7185.9,7124.26,6867.93,6466.64},
{7041.93,7066.04,7121.88,7185.06,7232.75,7226.47,7261.18,7243.58,7057,6753.06},
{7036,7052.97,7101.73,7180.63,7235.16,7269.31,7322.45,7333.18,7224.05,6980.21},
{7025.32,7035.4,7081.18,7164.75,7225.01,7294.92,7354.47,7391.9,7352.35,7172.36},
{7015.75,7025.41,7062.97,7150.88,7207.05,7298.16,7364.29,7425.26,7439.25,7347.4},
{7005.23,7013.6,7042.44,7125.77,7180.47,7292.44,7363.93,7447.99,7524.66,7479.74},
{6995.72,7016.61,7026.45,7106.14,7154.66,7275.21,7355.14,7450.69,7570.79,7587.36},
{6990.84,6990.99,7007.75,7078.07,7125.59,7247.77,7326.76,7438.9,7606.88,7680.67},
{6979.96,6977.98,6994.97,7053.51,7096.05,7219.43,7296.95,7415.57,7627.25,7751.46},
{6966.55,6969.74,6975.92,7007.94,7062.69,7148.28,7264.67,7390.13,7627.3,7799.28},
{6955.43,6957.7,6967.73,6961.1,7035.85,7110.7,7223.23,7354.83,7609.96,7824.3},
{6942.66,6942.68,6953.58,6969.7,7008.82,7077.73,7180.01,7308.72,7587.53,7841.62}};

	const double r4 = 0.004;
	const double r6 = 0.006075;
	const double r9 = 0.009;

	const Real burnup_avg_ = _burnup_avg[_qp];
	const Real density_ = _density[_qp];

	//ensuring the pellet radius falls within the range of numbers in the array for pellet_radius	
	if (_pellet_radius_property < r4)
	{
		mooseError("Pellet radius is less than 4 - value is outside of the array");
	}
	else if (_pellet_radius_property > r9)
	{
		mooseError("Pellet radius is greater than 9 - value is outside of the array");
	}
	//ensuring the values fit on the array
	if (burnup_avg_ > 960)
	{
		mooseError("Burnup average is greather than 960 - value is outside of the array");
	}
	if (burnup_avg_ < 0)
	{
		mooseError("burnup average is less than 0 - value is outside of the array and is non-physical");
	}
	if (_enrichment_property > 6)
	{
		mooseError("_enrichment is greather than 6 - value is outside of the array");
	}
	if (_enrichment_property < 0.71)
	{
		mooseError("_enrichment value is less than 0.71 - value is outside of the array");
	}

	double r4_val;
	double r6_val;
	double r9_val;

	double _beta;
	double _kappa;
	double _lambda;

	//interpolate between the r=4 and r=6 arrays if the pellet_radius is between those two values
	if (_pellet_radius_property < r6) //if the radius is between 4 and 6.075 mm
	{
		//finding the beta values on each of the two separate arrays
		r4_val = FissionHeatMaterial::interp2(r4_beta);
		r6_val = FissionHeatMaterial::interp2(r6_beta);
		//interpolating across the two arrays
		_beta = (r6_val - r4_val)/(r6 - r4) * (_pellet_radius_property - r4) + r4_val;

		//finding the kappa values on each of the two separate arrays
		r4_val = FissionHeatMaterial::interp2(r4_kappa);
		r6_val = FissionHeatMaterial::interp2(r6_kappa);
		//interpolating across the two arrays
		_kappa = (r6_val - r4_val)/(r6 - r4) * (_pellet_radius_property - r4) + r4_val;

		//finding the lambda values on each of the two separate arrays
		r4_val = FissionHeatMaterial::interp2(r4_lambda);
		r6_val = FissionHeatMaterial::interp2(r6_lambda);
		//interpolating across the two arrays
		_lambda = (r6_val - r4_val)/(r6 - r4) * (_pellet_radius_property - r4) + r4_val;
	}
	else if (_pellet_radius_property >= r6 && _pellet_radius_property < r9)
	{
		//finding the beta values on each of the two separate arrays
		r6_val = FissionHeatMaterial::interp2(r6_beta);
		r9_val = FissionHeatMaterial::interp2(r9_beta);
		//interpolating across the two arrays
		_beta = (r9_val - r6_val)/(r9 - r6) * (_pellet_radius_property - r6) + r6_val;

		//finding the kappa values on each of the two separate arrays
		r6_val = FissionHeatMaterial::interp2(r6_kappa);
		r9_val = FissionHeatMaterial::interp2(r9_kappa);
		//interpolating across the two arrays
		_kappa = (r9_val - r6_val)/(r9 - r6) * (_pellet_radius_property - r6) + r6_val;

		//finding the lambda values on each of the two separate arrays
		r6_val = FissionHeatMaterial::interp2(r6_lambda);
		r9_val = FissionHeatMaterial::interp2(r9_lambda);
		//interpolating across the two arrays
		_lambda = (r9_val - r6_val)/(r9 - r6) * (_pellet_radius_property - r6) + r6_val;
	}
	else if (_pellet_radius_property == r9)
	{
		_beta = FissionHeatMaterial::interp2(r9_beta);
		_kappa = FissionHeatMaterial::interp2(r9_kappa);
		_lambda = FissionHeatMaterial::interp2(r9_lambda);
	}

	Real r = 0;

	const Real fmag = (_linear_power->value(_t, _q_point[_qp])/(2 * pi)) * (1 / ((_pellet_radius_property / _kappa) * bessel1(_kappa * _pellet_radius_property) + _beta * ((_pellet_radius_property / _lambda) + (1 / std::pow(_lambda, 2)) * (std::exp(-_pellet_radius_property * _lambda) - 1)))); //eq: 5.23 from Prudil
	const Real vol_factor = density_ / _initial_density; //The change in volumetric heat generation due to thermal expansion and volume change
	
//If the geometry is 3D, use x,y,z coordinates - and compute r using the pythagorean theorem
	if (_is_3D)
	{
		Real x = _q_point[_qp](0);
		Real y = _q_point[_qp](1);
		r = std::sqrt(x * x + y * y);
	}
	//If the geometry is 2D, simply use and r-z coordinate plane where "x" in Trelis is considered "r" by MOOSE
	else if (!_is_3D)
	{
		r = _q_point[_qp](0);
	}
	Real Q = 0;
//eq: 5.20 from Prudil (with the middle term set to zero),
	if (_model_plate_fuel == false)
	{ 
		Q = vol_factor * fmag * (bessel0(_kappa * r) + _beta * std::exp(_lambda * (r - _pellet_radius_property))); 
	}
	else if (_model_plate_fuel == true)
	{
		Q = _linear_power->value(_t, _q_point[_qp]) / _area_fuel;
	}
	return Q; //the volumetric heat generation
}

//Defining the initial conditions
void
FissionHeatMaterial::initQpStatefulProperties()
{
	_q_fission[_qp] = _initial_qfission;
	_ratio[_qp] = _ratioProperty;
	_enrichment[_qp] = _enrichment_property;
	_pellet_radius[_qp] = _pellet_radius_property;
}

void
FissionHeatMaterial::computeQpProperties()
{
	if (_model_Qfission)
	{
		for(unsigned int qp(0); qp < _qrule->n_points(); ++qp)
  	{  
			Real Q(0);
			Q = FissionHeatMaterial::Q_fission();
			_q_fission[_qp] = Q;
  	}
	}
}
