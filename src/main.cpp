#include <iostream>
#include <vector>
#include <cmath>

struct State
{
	double u;
	double v;
	double P;
	double T;

	State operator+(const State& rhs) const
	{
		return {rhs.u+u, rhs.v+v, rhs.P+P, rhs.T+T};
	}
	State operator-(const State& rhs) const
	{	
		return {u-rhs.u, v-rhs.v, P-rhs.P, T-rhs.T};
	}
	State operator/(double rhs) const
	{
		return {u/rhs, v/rhs, P/rhs, T/rhs};
	}

	State(double u, double v, double P, double T):u(u), v(v), P(P), T(T){};
	State():u(0.0), v(0.0), P(0.0), T(0.0){};
};

void CalcExtrapolationsBC(std::vector<State>& field, size_t NPOINTSX, size_t NPOINTSY);

double CalcDensity(double P, double R, double T);
double CalcViscSutherland(double mu0, double T0, double T);

int main()
{
	const double machRef = 4.0;
	const double plateLength = 0.00001; // 0,01 mm
	const double aSeaLevel = 340.28;
	const double velRef = machRef * aSeaLevel;
	const double pressureSeaLevel = 101325;
	const double temperatureSeaLevel = 288.16;

	const double gama = 1.4; // Cp/Cv
	const double prandtlNumber = 0.71;
	const double muRef = 1.7894e-5;
	const double tempRef = temperatureSeaLevel;
	const double R = 287.0;

	const double Cv = R / (gama - 1);
	const double Cp = gama * Cv;

	const double RHOref = pressureSeaLevel/(R*tempRef);

	const double KFudgeFactor = 0.6;

	const size_t NPOINTSX = 70;
	const size_t NPOINTSY = 70;
	const size_t NITERATIONS = 100;

	const double LHORI = plateLength;
	const double LVERT = 25.0 * plateLength/sqrt(RHOref*velRef*plateLength/muRef);

	const double dx = plateLength / (NPOINTSX -1);
	const double dy = LVERT / (NPOINTSY-1);

	std::vector<State> U;
	std::vector<State> Up1;	
	std::vector<State> Ub;
	std::vector<State> dUdtb;

	U.resize(NPOINTSY* NPOINTSX);
	
	// SET INITIAL AND BOUNDARY CONDITIONS

	State defState = {0.0, 0.0, pressureSeaLevel, tempRef};

	for(size_t i = 0; i < NPOINTSY; i++)
		for(size_t j = 0; j < NPOINTSX; j++)
			U[j + (i*NPOINTSX)] = defState;

	const double PINF = pressureSeaLevel;
	const double TINF = temperatureSeaLevel;
	const double TWALL =  tempRef;
	const double UINF = velRef;

	// CASE 1 (LEADING EDGE)
	U[(NPOINTSY - 1) + (0*NPOINTSX)].P = PINF;
	U[(NPOINTSY - 1) + (0*NPOINTSX)].T = TINF;
	U[(NPOINTSY - 1) + (0*NPOINTSX)].u = 0.0;
	U[(NPOINTSY - 1) + (0*NPOINTSX)].v = 0.0;

	// CASE 2 (INFLOW/LEFT.TOP)

	for(size_t i = 0; i < NPOINTSY; i++)
	{	
		U[0 + (i*NPOINTSX)].P = PINF;
		U[0 + (i*NPOINTSX)].T = TINF;
		U[0 + (i*NPOINTSX)].u = UINF;
		U[0 + (i*NPOINTSX)].v = 0.0;
	}

	for(size_t j = 0; j < NPOINTSX; j++)
	{
		U[j + (0*NPOINTSX)].P = PINF;
		U[j + (0*NPOINTSX)].T = TINF;
		U[j + (0*NPOINTSX)].u = UINF;
		U[j + (0*NPOINTSX)].v = 0.0;
	}

	// CASE 3 (SURFACE EXCEPT LEADING EDGE)

	for(size_t j = 1; j < NPOINTSX; j++)
	{
		U[j + ((NPOINTSY-1)*NPOINTSX)].P = PINF;
		U[j + ((NPOINTSY-1)*NPOINTSX)].T = TWALL;
		U[j + ((NPOINTSY-1)*NPOINTSX)].u = 0.0;
		U[j + ((NPOINTSY-1)*NPOINTSX)].v = 0.0;
	}

	CalcExtrapolationsBC(U, NPOINTSX, NPOINTSY);

	Up1 = U;
	Ub = U;

	
	for(size_t it = 0; it < NITERATIONS; it++)
	{
		std::cout << std::endl;
		std::cout << "Iteration number: " << it << std::endl;

		// Calculate time step;

		double dt = 1000.0;

		for(size_t i = 0; i < NPOINTSY; i++)
		{
			for(size_t j = 0; j < NPOINTSX; j++)
			{
				const State& Uij = U[j + (i*NPOINTSX)];

				const double aij = sqrt(gama * R * Uij.T);
				const double mu = CalcViscSutherland(muRef, tempRef, Uij.T);

				const double vlij =  (4.0 / 3.0) * mu * gama * (mu /prandtlNumber)/CalcDensity(Uij.P, R, Uij.T);
				const double suminverse = (1.0/(dx*dx)) + (1.0/(dy*dy));

				const double deltatcfl = KFudgeFactor/((abs(Uij.u)/dx) + (abs(Uij.v)/dy) + (aij*sqrt(suminverse)) 
				+ (2.0 * vlij*suminverse));

				if(deltatcfl < dt)
				{
					dt = deltatcfl;
				}
			}
		}

		std::cout << "dt: " << dt << std::endl;

		
		

	}

	
	return 0;
}

double CalcDensity(double P, double R, double T)
{
	return P/(R*T);
}
double CalcViscSutherland(double mu0, double T0, double T)
{
	return mu0*pow(T/T0, 1.5)*((T0+110.0)/(T+110.0));
}

void CalcExtrapolationsBC(std::vector<State>& field, size_t NPOINTSX, size_t NPOINTSY)
{
	// CASE 3 (surface except leadign edge)

	for(size_t j = 1; j < NPOINTSX; j++)
	{
		const double pextr = (2.0 * field[j + ((NPOINTSY-2)*NPOINTSX)].P) - 
		(field[j + ((NPOINTSY-3)*NPOINTSX)].P);

		field[j + ((NPOINTSY-1)*NPOINTSX)].P = pextr;
	}

	// CASE 4 (outflow except surface and inflow)

	for(size_t i = 1; i < NPOINTSY - 1; i++)
	{
		const double pextr = (2.0 * field[NPOINTSX - 2 + (i*NPOINTSX)].P) - 
		(field[NPOINTSX - 3 + (i*NPOINTSX)].P);
		const double uextr = (2.0 * field[NPOINTSX - 2 + (i*NPOINTSX)].u) - 
		(field[NPOINTSX - 3 + (i*NPOINTSX)].u);
		const double vextr = (2.0 * field[NPOINTSX - 2 + (i*NPOINTSX)].v) - 
		(field[NPOINTSX - 3 + (i*NPOINTSX)].v);
		const double textr = (2.0 * field[NPOINTSX - 2 + (i*NPOINTSX)].T) - 
		(field[NPOINTSX - 3 + (i*NPOINTSX)].T);

		field[NPOINTSX - 1 + (i*NPOINTSX)].P = pextr;
		field[NPOINTSX - 1 + (i*NPOINTSX)].u = uextr;
		field[NPOINTSX - 1 + (i*NPOINTSX)].v = vextr;
		field[NPOINTSX - 1 + (i*NPOINTSX)].T = textr;
	}	
}