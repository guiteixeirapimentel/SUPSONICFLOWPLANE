#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

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
	State operator*(double rhs) const
	{
		return {u*rhs, v*rhs, P*rhs, T*rhs};
	}
	State operator-() const
	{
		return {-u, -v, -P, -T};
	}

	State(double u, double v, double P, double T):u(u), v(v), P(P), T(T){};
	State():u(0.0), v(0.0), P(0.0), T(0.0){};
};

State CalcState(const State& U);
State CalcU(const State& estado);

void CalcExtrapolationsBC(std::vector<State>& field, size_t NPOINTSX, size_t NPOINTSY);

double CalcDensity(double P, double R, double T);
double CalcViscSutherland(double mu0, double T0, double T);
double CalcEt(double rho, double Cv, double T, double V);
double Calck(double mu);

// i j in "array notation" (change in row and collumn respectively)
double GetE1Predictor(const std::vector<State>& field, size_t i, size_t j, size_t NPOINTSX, size_t NPOINTSY);
// i j in "array notation" (change in row and collumn respectively)
double GetE2Predictor(const std::vector<State>& field, size_t i, size_t j, size_t NPOINTSX, size_t NPOINTSY);
// i j in "array notation" (change in row and collumn respectively)
double GetE3Predictor(const std::vector<State>& field, size_t i, size_t j, size_t NPOINTSX, size_t NPOINTSY);
// i j in "array notation" (change in row and collumn respectively)
double GetE4Predictor(const std::vector<State>& field, size_t i, size_t j, size_t NPOINTSX, size_t NPOINTSY);

// i j in "array notation" (change in row and collumn respectively)
double GetF1Predictor(const std::vector<State>& field, size_t i, size_t j, size_t NPOINTSX, size_t NPOINTSY);
// i j in "array notation" (change in row and collumn respectively)
double GetF2Predictor(const std::vector<State>& field, size_t i, size_t j, size_t NPOINTSX, size_t NPOINTSY);
// i j in "array notation" (change in row and collumn respectively)
double GetF3Predictor(const std::vector<State>& field, size_t i, size_t j, size_t NPOINTSX, size_t NPOINTSY);
// i j in "array notation" (change in row and collumn respectively)
double GetF4Predictor(const std::vector<State>& field, size_t i, size_t j, size_t NPOINTSX, size_t NPOINTSY);

// i j in "array notation" (change in row and collumn respectively)
double CalcTauXXPredictor(const std::vector<State>& field, size_t i, size_t j, size_t NPOINTSX, size_t NPOINTSY);
// i j in "array notation" (change in row and collumn respectively)
double CalcTauYYPredictor(const std::vector<State>& field, size_t i, size_t j, size_t NPOINTSX, size_t NPOINTSY);
// i j in "array notation" (change in row and collumn respectively)
double CalcTauXYPredictorE(const std::vector<State>& field, size_t i, size_t j, size_t NPOINTSX, size_t NPOINTSY);
// i j in "array notation" (change in row and collumn respectively)
double CalcTauXYPredictorF(const std::vector<State>& field, size_t i, size_t j, size_t NPOINTSX, size_t NPOINTSY);
// i j in "array notation" (change in row and collumn respectively)
double CalcQxPredictor(const std::vector<State>& field, size_t i, size_t j, size_t NPOINTSX, size_t NPOINTSY);
// i j in "array notation" (change in row and collumn respectively)
double CalcQyPredictor(const std::vector<State>& field, size_t i, size_t j, size_t NPOINTSX, size_t NPOINTSY);

// i j in "array notation" (change in row and collumn respectively)
double GetE1Corrector(const std::vector<State>& field, size_t i, size_t j, size_t NPOINTSX, size_t NPOINTSY);
// i j in "array notation" (change in row and collumn respectively)
double GetE2Corrector(const std::vector<State>& field, size_t i, size_t j, size_t NPOINTSX, size_t NPOINTSY);
// i j in "array notation" (change in row and collumn respectively)
double GetE3Corrector(const std::vector<State>& field, size_t i, size_t j, size_t NPOINTSX, size_t NPOINTSY);
// i j in "array notation" (change in row and collumn respectively)
double GetE4Corrector(const std::vector<State>& field, size_t i, size_t j, size_t NPOINTSX, size_t NPOINTSY);

// i j in "array notation" (change in row and collumn respectively)
double GetF1Corrector(const std::vector<State>& field, size_t i, size_t j, size_t NPOINTSX, size_t NPOINTSY);
// i j in "array notation" (change in row and collumn respectively)
double GetF2Corrector(const std::vector<State>& field, size_t i, size_t j, size_t NPOINTSX, size_t NPOINTSY);
// i j in "array notation" (change in row and collumn respectively)
double GetF3Corrector(const std::vector<State>& field, size_t i, size_t j, size_t NPOINTSX, size_t NPOINTSY);
// i j in "array notation" (change in row and collumn respectively)
double GetF4Corrector(const std::vector<State>& field, size_t i, size_t j, size_t NPOINTSX, size_t NPOINTSY);

// i j in "array notation" (change in row and collumn respectively)
double CalcTauXXCorrector(const std::vector<State>& field, size_t i, size_t j, size_t NPOINTSX, size_t NPOINTSY);
// i j in "array notation" (change in row and collumn respectively)
double CalcTauYYCorrector(const std::vector<State>& field, size_t i, size_t j, size_t NPOINTSX, size_t NPOINTSY);
// i j in "array notation" (change in row and collumn respectively)
double CalcTauXYCorrectorE(const std::vector<State>& field, size_t i, size_t j, size_t NPOINTSX, size_t NPOINTSY);
// i j in "array notation" (change in row and collumn respectively)
double CalcTauXYCorrectorF(const std::vector<State>& field, size_t i, size_t j, size_t NPOINTSX, size_t NPOINTSY);
// i j in "array notation" (change in row and collumn respectively)
double CalcQxCorrector(const std::vector<State>& field, size_t i, size_t j, size_t NPOINTSX, size_t NPOINTSY);
// i j in "array notation" (change in row and collumn respectively)
double CalcQyCorrector(const std::vector<State>& field, size_t i, size_t j, size_t NPOINTSX, size_t NPOINTSY);

constexpr double machRef = 4.0;
constexpr double plateLength = 0.00001; // 0,01 mm
constexpr double aSeaLevel = 340.28;
constexpr double velRef = machRef * aSeaLevel;
constexpr double pressureSeaLevel = 101325;
constexpr double temperatureSeaLevel = 288.16;
constexpr double gama = 1.4; // Cp/Cv
constexpr double prandtlNumber = 0.71;
constexpr double muRef = 1.7894e-5;
constexpr double tempRef = temperatureSeaLevel;
constexpr double R = 287.0;
constexpr double Cv = R / (gama - 1);
constexpr double Cp = gama * Cv;
constexpr double RHOref = pressureSeaLevel/(R*tempRef);
constexpr double KFudgeFactor = 0.6;
constexpr size_t NPOINTSX = 70;
constexpr size_t NPOINTSY = 70;
constexpr size_t NITERATIONS = 100;
constexpr double LHORI = plateLength;
const double LVERT = 25.0 * plateLength/sqrt(RHOref*velRef*plateLength/muRef);
constexpr double dx = plateLength / (NPOINTSX -1);
const double dy = LVERT / (NPOINTSY-1);

int main()
{
	std::vector<State> U;
	std::vector<State> estado;
	std::vector<State> estadoPredicted;
	std::vector<State> Upredicted;	
	std::vector<State> Ucorrected;
	std::vector<State> dUdtpredicted;
	std::vector<State> dUdtcorrected;

	estado.resize(NPOINTSY* NPOINTSX);
	
	// SET INITIAL AND BOUNDARY CONDITIONS

	State defState = {0.0, 0.0, pressureSeaLevel, tempRef};

	for(size_t i = 0; i < NPOINTSY; i++)
		for(size_t j = 0; j < NPOINTSX; j++)
			estado[j + (i*NPOINTSX)] = defState;

	const double PINF = pressureSeaLevel;
	const double TINF = temperatureSeaLevel;
	const double TWALL =  tempRef;
	const double UINF = velRef;

	// CASE 1 (LEADING EDGE)
	estado[(NPOINTSY - 1) + (0*NPOINTSX)].P = PINF;
	estado[(NPOINTSY - 1) + (0*NPOINTSX)].T = TINF;
	estado[(NPOINTSY - 1) + (0*NPOINTSX)].u = 0.0;
	estado[(NPOINTSY - 1) + (0*NPOINTSX)].v = 0.0;

	// CASE 2 (INFLOW/LEFT.TOP)

	for(size_t i = 0; i < NPOINTSY; i++)
	{	
		estado[0 + (i*NPOINTSX)].P = PINF;
		estado[0 + (i*NPOINTSX)].T = TINF;
		estado[0 + (i*NPOINTSX)].u = UINF;
		estado[0 + (i*NPOINTSX)].v = 0.0;
	}

	for(size_t j = 0; j < NPOINTSX; j++)
	{
		estado[j + (0*NPOINTSX)].P = PINF;
		estado[j + (0*NPOINTSX)].T = TINF;
		estado[j + (0*NPOINTSX)].u = UINF;
		estado[j + (0*NPOINTSX)].v = 0.0;
	}

	// CASE 3 (SURFACE EXCEPT LEADING EDGE)

	for(size_t j = 1; j < NPOINTSX; j++)
	{
		estado[j + ((NPOINTSY-1)*NPOINTSX)].P = PINF;
		estado[j + ((NPOINTSY-1)*NPOINTSX)].T = TWALL;
		estado[j + ((NPOINTSY-1)*NPOINTSX)].u = 0.0;
		estado[j + ((NPOINTSY-1)*NPOINTSX)].v = 0.0;
	}

	CalcExtrapolationsBC(estado, NPOINTSX, NPOINTSY);

	estadoPredicted = estado;

	U.resize(NPOINTSX*NPOINTSY);

	for(size_t i =0; i < NPOINTSY; i++)
	{
		for(size_t j = 0; j < NPOINTSX; j++)
		{
			U[j + (i*NPOINTSX)] = CalcU(estado[j + (i*NPOINTSX)]);
		}
	}

	Upredicted = U;
	Ucorrected = U;

	dUdtpredicted = U;
	dUdtcorrected = U;

	double totalTime = 0.0;
	
	for(size_t it = 0; it < NITERATIONS; it++)
	{
		std::cout << std::endl;
		std::cout << "Iteration number: " << it << std::endl;

		// Calculate time step;

		double dt = 1000.0;

		double maxVlij = 0.0;

		for(size_t i = 0; i < NPOINTSY; i++)
		{
			for(size_t j = 0; j < NPOINTSX; j++)
			{
				const State& estadoij = estado[j + (i*NPOINTSX)];

				const double mu = CalcViscSutherland(muRef, tempRef, estadoij.T);

				const double vlij =  (4.0 / 3.0) * mu * gama * (mu /prandtlNumber)/CalcDensity(estadoij.P, R, estadoij.T);

				if(vlij > maxVlij)
					maxVlij = vlij;
			}
		}

		for(size_t i = 0; i < NPOINTSY; i++)
		{
			for(size_t j = 0; j < NPOINTSX; j++)
			{
				const State& estadoij = estado[j + (i*NPOINTSX)];

				const double aij = sqrt(gama * R * estadoij.T);

				const double suminverse = (1.0/(dx*dx)) + (1.0/(dy*dy));

				const double deltatcfl = KFudgeFactor/((fabs(estadoij.u)/dx) + (fabs(estadoij.v)/dy) + (aij*sqrt(suminverse)) 
				+ (2.0 * maxVlij*suminverse));

				if(deltatcfl < dt)
				{
					dt = deltatcfl;
				}
			}
		}

		std::cout << "dt: " << dt << std::endl;

		// CALCULATE ALL PREDICTED DERIVATIVES

		// Plus in book -> minus in the array.
		for(size_t i = 1; i < NPOINTSY - 1; i++)
		{
			for(size_t j = 1; j < NPOINTSX - 1; j++)
			{				
				const double E1ipj = GetE1Predictor(estado, i, j+1, NPOINTSX, NPOINTSY);
				const double E2ipj = GetE2Predictor(estado, i, j+1, NPOINTSX, NPOINTSY);
				const double E3ipj = GetE3Predictor(estado, i, j+1, NPOINTSX, NPOINTSY);
				const double E4ipj = GetE4Predictor(estado, i, j+1, NPOINTSX, NPOINTSY);

				const double E1ij = GetE1Predictor(estado, i, j, NPOINTSX, NPOINTSY);
				const double E2ij = GetE2Predictor(estado, i, j, NPOINTSX, NPOINTSY);
				const double E3ij = GetE3Predictor(estado, i, j, NPOINTSX, NPOINTSY);
				const double E4ij = GetE4Predictor(estado, i, j, NPOINTSX, NPOINTSY);
				
				const double F1ijp = GetF1Predictor(estado, i-1, j, NPOINTSX, NPOINTSY);
				const double F2ijp = GetF2Predictor(estado, i-1, j, NPOINTSX, NPOINTSY);
				const double F3ijp = GetF3Predictor(estado, i-1, j, NPOINTSX, NPOINTSY);
				const double F4ijp = GetF4Predictor(estado, i-1, j, NPOINTSX, NPOINTSY);

				const double F1ij = GetF1Predictor(estado, i, j, NPOINTSX, NPOINTSY);
				const double F2ij = GetF2Predictor(estado, i, j, NPOINTSX, NPOINTSY);
				const double F3ij = GetF3Predictor(estado, i, j, NPOINTSX, NPOINTSY);
				const double F4ij = GetF4Predictor(estado, i, j, NPOINTSX, NPOINTSY);

				const State Eipj = {E1ipj, E2ipj, E3ipj, E4ipj};
				const State Eij = {E1ij, E2ij, E3ij, E4ij};

				const State Fijp = {F1ijp, F2ijp, F3ijp, F4ijp};
				const State Fij = {F1ij, F2ij, F3ij, F4ij};

				dUdtpredicted[j + (i * NPOINTSX)] = -( ((Eipj-Eij)/dx) + ((Fijp - Fij)/dy) );				
			}
		}

		if(it == 58)
		{
			float x = 0;
		}
		// CALCULATE ALL THE PREDICTED VALUES

		for(size_t i = 1; i < NPOINTSY - 1; i++)
		{
			for(size_t j = 1; j < NPOINTSX - 1; j++)
			{
				if(j == 1)
				{
					float x = 0;
				}
				Upredicted[j + (i*NPOINTSX)] = U[j + (i*NPOINTSX)] + 
				(dUdtpredicted[j + (i * NPOINTSX)] * dt);

				estadoPredicted[j + (i*NPOINTSX)] = CalcState(Upredicted[j + (i*NPOINTSX)]);

				if(estadoPredicted[j+(i*NPOINTSX)].T < 0.0)
				{
					std::cout << "Erro i: " << i << " j: " << j << " on predicted values" << std::endl;
				}
			}
		}

		CalcExtrapolationsBC(estadoPredicted, NPOINTSX, NPOINTSY);
		
		// CALCULATE ALL THE CORRECTED DERIVATIVES

		for(size_t i = 1; i < NPOINTSY - 1; i++)
		{
			for(size_t j = 1; j < NPOINTSX - 1; j++)
			{				
				const double E1imj = GetE1Corrector(estadoPredicted, i, j-1, NPOINTSX, NPOINTSY);
				const double E2imj = GetE2Corrector(estadoPredicted, i, j-1, NPOINTSX, NPOINTSY);
				const double E3imj = GetE3Corrector(estadoPredicted, i, j-1, NPOINTSX, NPOINTSY);
				const double E4imj = GetE4Corrector(estadoPredicted, i, j-1, NPOINTSX, NPOINTSY);

				const double E1ij = GetE1Corrector(estadoPredicted, i, j, NPOINTSX, NPOINTSY);
				const double E2ij = GetE2Corrector(estadoPredicted, i, j, NPOINTSX, NPOINTSY);
				const double E3ij = GetE3Corrector(estadoPredicted, i, j, NPOINTSX, NPOINTSY);
				const double E4ij = GetE4Corrector(estadoPredicted, i, j, NPOINTSX, NPOINTSY);
				
				const double F1ijm = GetF1Corrector(estadoPredicted, i+1, j, NPOINTSX, NPOINTSY);
				const double F2ijm = GetF2Corrector(estadoPredicted, i+1, j, NPOINTSX, NPOINTSY);
				const double F3ijm = GetF3Corrector(estadoPredicted, i+1, j, NPOINTSX, NPOINTSY);
				const double F4ijm = GetF4Corrector(estadoPredicted, i+1, j, NPOINTSX, NPOINTSY);

				const double F1ij = GetF1Corrector(estadoPredicted, i, j, NPOINTSX, NPOINTSY);
				const double F2ij = GetF2Corrector(estadoPredicted, i, j, NPOINTSX, NPOINTSY);
				const double F3ij = GetF3Corrector(estadoPredicted, i, j, NPOINTSX, NPOINTSY);
				const double F4ij = GetF4Corrector(estadoPredicted, i, j, NPOINTSX, NPOINTSY);

				const State Eimj = {E1ij, E2imj, E3imj, E4imj};
				const State Eij = {E1ij, E2ij, E3ij, E4ij};

				const State Fijm = {F1ijm, F2ijm, F3ijm, F4ijm};
				const State Fij = {F1ij, F2ij, F3ij, F4ij};

				dUdtcorrected[j + (i * NPOINTSX)] = -( ((Eij-Eimj)/dx) + ((Fij - Fijm)/dy) );				
			}
		}

		// CALCULATE ALL THE CORRECTED VALUES
		
		for(size_t i = 1; i < NPOINTSY - 1; i++)
		{
			for(size_t j = 1; j < NPOINTSX - 1; j++)
			{
				Ucorrected[j + (i*NPOINTSX)] = U[j + (i*NPOINTSX)] + 
				((dUdtpredicted[j + (i * NPOINTSX)] + dUdtcorrected[j + (i * NPOINTSX)])*(dt/2.0));

				estado[j + (i*NPOINTSX)] = CalcState(Ucorrected[j + (i*NPOINTSX)]);

				if(estado[j+(i*NPOINTSX)].T < 0.0)
				{
					std::cout << "Erro i: " << i << " j: " << j << " on corrected values" << std::endl;
				}

			}
		}

		U = Ucorrected;

		CalcExtrapolationsBC(estado, NPOINTSX, NPOINTSY);

		totalTime += dt;
		std::cout << "Total time passed: " << totalTime << std::endl;
	}

	std::ofstream arqu("u.csv");

	for(size_t i = 0; i < NPOINTSY; i++)
	{
		for(size_t j = 0; j < NPOINTSX; j++)
		{
			arqu << estado[j + (i * NPOINTSX)].u << ";";

		}
		arqu << "\n";
	}

	arqu.close();

	std::ofstream arqv("v.csv");

	for(size_t i = 0; i < NPOINTSY; i++)
	{
		for(size_t j = 0; j < NPOINTSX; j++)
		{
			arqv << estado[j + (i * NPOINTSX)].v << ";";

		}
		arqv << "\n";
	}

	arqv.close();

	std::ofstream arqp("P.csv");

	for(size_t i = 0; i < NPOINTSY; i++)
	{
		for(size_t j = 0; j < NPOINTSX; j++)
		{
			arqp << estado[j + (i * NPOINTSX)].P << ";";

		}
		arqp << "\n";
	}

	arqp.close();

	std::ofstream arqT("T.csv");

	for(size_t i = 0; i < NPOINTSY; i++)
	{
		for(size_t j = 0; j < NPOINTSX; j++)
		{
			arqT << estado[j + (i * NPOINTSX)].T << ";";

		}
		arqT << "\n";
	}

	arqT.close();
	
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

double Calck(double mu)
{
	return mu*Cp/prandtlNumber;
}

double CalcEt(double rho, double Cv, double T, double Vquad)
{
	return rho*((Cv * T) + (Vquad/2.0));
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

// i j in "array notation" (change in row and collumn respectively)
double GetE1Predictor(const std::vector<State>& field, size_t i, size_t j, size_t NPOINTSX, size_t NPOINTSY)
{
	const State& fieldij = field[j + (i * NPOINTSX)];
	const double res = fieldij.u * CalcDensity(fieldij.P, R, fieldij.T);

	return res;
}

// i j in "array notation" (change in row and collumn respectively)
double GetE2Predictor(const std::vector<State>& field, size_t i, size_t j, size_t NPOINTSX, size_t NPOINTSY)
{
	const State& fieldij = field[j + (i * NPOINTSX)];
	double res = fieldij.u * fieldij.u * CalcDensity(fieldij.P, R, fieldij.T);

	res += fieldij.P;
	res -= CalcTauXXPredictor(field, i, j, NPOINTSX, NPOINTSY);

	return res;
}
// i j in "array notation" (change in row and collumn respectively)
double GetE3Predictor(const std::vector<State>& field, size_t i, size_t j, size_t NPOINTSX, size_t NPOINTSY)
{
	const State& fieldij = field[j + (i * NPOINTSX)];
	double res = fieldij.v * fieldij.u * CalcDensity(fieldij.P, R, fieldij.T);

	res -= CalcTauXYPredictorE(field, i, j, NPOINTSX, NPOINTSY);

	return res;
}
// i j in "array notation" (change in row and collumn respectively)
double GetE4Predictor(const std::vector<State>& field, size_t i, size_t j, size_t NPOINTSX, size_t NPOINTSY)
{
	const State& fieldij = field[j + (i * NPOINTSX)];
	const double Vquad = (fieldij.u*fieldij.u) + (fieldij.v*fieldij.v);
	const double rho = CalcDensity(fieldij.P, R, fieldij.T);
	double res = fieldij.u * CalcEt(rho, Cv, fieldij.T, Vquad);
	res += fieldij.P*fieldij.u;
	res -= fieldij.u * CalcTauXXPredictor(field, i, j, NPOINTSX, NPOINTSY);
	res -= fieldij.v * CalcTauXYPredictorE(field, i, j, NPOINTSX, NPOINTSY);
	res += CalcQxPredictor(field, i, j, NPOINTSX, NPOINTSY);

	return res;
}

// i j in "array notation" (change in row and collumn respectively)
double GetF1Predictor(const std::vector<State>& field, size_t i, size_t j, size_t NPOINTSX, size_t NPOINTSY)
{
	const State& fieldij = field[j + (i * NPOINTSX)];
	const double res = fieldij.v * CalcDensity(fieldij.P, R, fieldij.T);

	return res;
}
// i j in "array notation" (change in row and collumn respectively)
double GetF2Predictor(const std::vector<State>& field, size_t i, size_t j, size_t NPOINTSX, size_t NPOINTSY)
{	
	const State& fieldij = field[j + (i * NPOINTSX)];
	double res = fieldij.u*fieldij.v * CalcDensity(fieldij.P, R, fieldij.T);

	res -= CalcTauXYPredictorF(field, i, j, NPOINTSX, NPOINTSY);

	return res;
}
// i j in "array notation" (change in row and collumn respectively)
double GetF3Predictor(const std::vector<State>& field, size_t i, size_t j, size_t NPOINTSX, size_t NPOINTSY)
{		
	const State& fieldij = field[j + (i * NPOINTSX)];
	double res = fieldij.v*fieldij.v * CalcDensity(fieldij.P, R, fieldij.T);
	res += fieldij.P;
	
	res -= CalcTauYYPredictor(field, i, j, NPOINTSX, NPOINTSY);

	return res;
}
// i j in "array notation" (change in row and collumn respectively)
double GetF4Predictor(const std::vector<State>& field, size_t i, size_t j, size_t NPOINTSX, size_t NPOINTSY)
{
	const State& fieldij = field[j + (i * NPOINTSX)];
	const double Vquad = (fieldij.u*fieldij.u) + (fieldij.v*fieldij.v);

	double res = fieldij.v * CalcEt(CalcDensity(fieldij.P, R, fieldij.T), Cv, fieldij.T, Vquad);
	res += fieldij.P*fieldij.v;
	res -= fieldij.u * CalcTauXYPredictorF(field, i, j, NPOINTSX, NPOINTSY);
	res -= fieldij.v * CalcTauYYPredictor(field, i, j, NPOINTSX, NPOINTSY);
	res += CalcQyPredictor(field, i, j, NPOINTSX, NPOINTSY);

	return res;
}

// i j in "array notation" (change in row and collumn respectively)
double CalcTauXXPredictor(const std::vector<State>& field, size_t i, size_t j, size_t NPOINTSX, size_t NPOINTSY)
{
	const State& fieldij = field[j + (i * NPOINTSX)];
	const State& fieldimj = field[j-1 + (i * NPOINTSX)];
	const State& fieldijp = field[j + ((i-1) * NPOINTSX)];
	const State& fieldijm = field[j + ((i+1) * NPOINTSX)];
	
	const double mu = CalcViscSutherland(muRef, tempRef, fieldij.T);

	const double res = (2.0/3.0) * mu * ( 
		((2.0/dx)*(fieldij.u - fieldimj.u)) - (((fieldijp.v - fieldijm.v)/(2.0*dy))) 
		);

	return res;
}
// i j in "array notation" (change in row and collumn respectively)
double CalcTauYYPredictor(const std::vector<State>& field, size_t i, size_t j, size_t NPOINTSX, size_t NPOINTSY)
{
	const State& fieldij = field[j + (i * NPOINTSX)];
	const State& fieldijm = field[j + ((i+1) * NPOINTSX)];
	const State& fieldipj = field[j+1 + (i * NPOINTSX)];
	const State& fieldimj = field[j-1 + (i * NPOINTSX)];
	
	const double mu = CalcViscSutherland(muRef, tempRef, fieldij.T);	

	const double res = (2.0/3.0) * mu * ( 
		((2.0/dy)*(fieldij.v - fieldijm.v)) - (((fieldipj.u - fieldimj.u)/(2.0*dx)))
		);

	return res;
}
// i j in "array notation" (change in row and collumn respectively)
double CalcTauXYPredictorE(const std::vector<State>& field, size_t i, size_t j, size_t NPOINTSX, size_t NPOINTSY)
{
	const State& fieldij = field[j + (i * NPOINTSX)];
	const State& fieldijp = field[j + ((i-1) * NPOINTSX)];
	const State& fieldijm = field[j + ((i+1) * NPOINTSX)];
	const State& fieldimj = field[j-1 + (i * NPOINTSX)];
	
	const double mu = CalcViscSutherland(muRef, tempRef, fieldij.T);

	const double res = mu * ( 
		((fieldij.v - fieldimj.v)/dx) + (((fieldijp.u - fieldijm.u)/(2.0*dy))) 
		);

	return res;
}
// i j in "array notation" (change in row and collumn respectively)
double CalcTauXYPredictorF(const std::vector<State>& field, size_t i, size_t j, size_t NPOINTSX, size_t NPOINTSY)
{
	const State& fieldij = field[j + (i * NPOINTSX)];
	const State& fieldijm = field[j + ((i+1) * NPOINTSX)];
	const State& fieldipj = field[j+1 + (i * NPOINTSX)];
	const State& fieldimj = field[j-1 + (i * NPOINTSX)];
	
	const double mu = CalcViscSutherland(muRef, tempRef, fieldij.T);

	const double res = mu * ( 
		((fieldij.u - fieldijm.u)/dy) + (((fieldipj.v - fieldimj.v)/(2.0*dx))) 
		);

	return res;
}
// i j in "array notation" (change in row and collumn respectively)
double CalcQxPredictor(const std::vector<State>& field, size_t i, size_t j, size_t NPOINTSX, size_t NPOINTSY)
{
	const State& fieldij = field[j + (i * NPOINTSX)];
	const State& fieldimj = field[j-1 + (i * NPOINTSX)];
	
	const double mu = CalcViscSutherland(muRef, tempRef, fieldij.T);
	const double k = Calck(mu);
	
	const double res = -k * (fieldij.T - fieldimj.T)/dx;

	return res;
}
// i j in "array notation" (change in row and collumn respectively)
double CalcQyPredictor(const std::vector<State>& field, size_t i, size_t j, size_t NPOINTSX, size_t NPOINTSY)
{
	const State& fieldij = field[j + (i * NPOINTSX)];
	const State& fieldijm = field[j + ((i+1) * NPOINTSX)];
	
	const double mu = CalcViscSutherland(muRef, tempRef, fieldij.T);
	const double k = Calck(mu);
	
	const double res = -k * (fieldij.T - fieldijm.T)/dy;

	return res;
}

// i j in "array notation" (change in row and collumn respectively)
double GetE1Corrector(const std::vector<State>& field, size_t i, size_t j, size_t NPOINTSX, size_t NPOINTSY)
{
	const State& fieldij = field[j + (i * NPOINTSX)];
	const double res = fieldij.u * CalcDensity(fieldij.P, R, fieldij.T);

	return res;
}

// i j in "array notation" (change in row and collumn respectively)
double GetE2Corrector(const std::vector<State>& field, size_t i, size_t j, size_t NPOINTSX, size_t NPOINTSY)
{
	const State& fieldij = field[j + (i * NPOINTSX)];
	double res = fieldij.u * fieldij.u * CalcDensity(fieldij.P, R, fieldij.T);

	res += fieldij.P;
	res -= CalcTauXXCorrector(field, i, j, NPOINTSX, NPOINTSY);

	return res;
}
// i j in "array notation" (change in row and collumn respectively)
double GetE3Corrector(const std::vector<State>& field, size_t i, size_t j, size_t NPOINTSX, size_t NPOINTSY)
{
	const State& fieldij = field[j + (i * NPOINTSX)];
	double res = fieldij.v * fieldij.u * CalcDensity(fieldij.P, R, fieldij.T);

	res -= CalcTauXYCorrectorE(field, i, j, NPOINTSX, NPOINTSY);

	return res;
}
// i j in "array notation" (change in row and collumn respectively)
double GetE4Corrector(const std::vector<State>& field, size_t i, size_t j, size_t NPOINTSX, size_t NPOINTSY)
{
	const State& fieldij = field[j + (i * NPOINTSX)];
	const double Vquad = (fieldij.u*fieldij.u) + (fieldij.v*fieldij.v);
	double res = fieldij.u * CalcEt(CalcDensity(fieldij.P, R, fieldij.T), Cv, fieldij.T, Vquad);
	res += fieldij.P*fieldij.u;
	res -= fieldij.u * CalcTauXXCorrector(field, i, j, NPOINTSX, NPOINTSY);
	res -= fieldij.v * CalcTauXYCorrectorE(field, i, j, NPOINTSX, NPOINTSY);
	res += CalcQxCorrector(field, i, j, NPOINTSX, NPOINTSY);

	return res;
}

// i j in "array notation" (change in row and collumn respectively)
double GetF1Corrector(const std::vector<State>& field, size_t i, size_t j, size_t NPOINTSX, size_t NPOINTSY)
{
	const State& fieldij = field[j + (i * NPOINTSX)];
	const double res = fieldij.v * CalcDensity(fieldij.P, R, fieldij.T);

	return res;
}
// i j in "array notation" (change in row and collumn respectively)
double GetF2Corrector(const std::vector<State>& field, size_t i, size_t j, size_t NPOINTSX, size_t NPOINTSY)
{	
	const State& fieldij = field[j + (i * NPOINTSX)];
	double res = fieldij.u*fieldij.v * CalcDensity(fieldij.P, R, fieldij.T);

	res -= CalcTauXYCorrectorF(field, i, j, NPOINTSX, NPOINTSY);

	return res;
}
// i j in "array notation" (change in row and collumn respectively)
double GetF3Corrector(const std::vector<State>& field, size_t i, size_t j, size_t NPOINTSX, size_t NPOINTSY)
{		
	const State& fieldij = field[j + (i * NPOINTSX)];
	double res = fieldij.v*fieldij.v * CalcDensity(fieldij.P, R, fieldij.T);
	res += fieldij.P;
	
	res -= CalcTauYYCorrector(field, i, j, NPOINTSX, NPOINTSY);

	return res;
}
// i j in "array notation" (change in row and collumn respectively)
double GetF4Corrector(const std::vector<State>& field, size_t i, size_t j, size_t NPOINTSX, size_t NPOINTSY)
{
	const State& fieldij = field[j + (i * NPOINTSX)];
	const double Vquad = (fieldij.u*fieldij.u) + (fieldij.v*fieldij.v);

	double res = fieldij.v * CalcEt(CalcDensity(fieldij.P, R, fieldij.T), Cv, fieldij.T, Vquad);
	res += fieldij.P*fieldij.v;
	res -= fieldij.u * CalcTauXYCorrectorF(field, i, j, NPOINTSX, NPOINTSY);
	res -= fieldij.v * CalcTauYYCorrector(field, i, j, NPOINTSX, NPOINTSY);
	res += CalcQyCorrector(field, i, j, NPOINTSX, NPOINTSY);

	return res;
}

// i j in "array notation" (change in row and collumn respectively)
double CalcTauXXCorrector(const std::vector<State>& field, size_t i, size_t j, size_t NPOINTSX, size_t NPOINTSY)
{
	const State& fieldij = field[j + (i * NPOINTSX)];
	const State& fieldipj = field[j+1 + (i * NPOINTSX)];
	const State& fieldijp = field[j + ((i-1) * NPOINTSX)];
	const State& fieldijm = field[j + ((i+1) * NPOINTSX)];
	
	const double mu = CalcViscSutherland(muRef, tempRef, fieldij.T);

	const double res = (2.0/3.0) * mu * ( 
		((2.0/dx)*(fieldipj.u - fieldij.u)) - (((fieldijp.v - fieldijm.v)/(2.0*dy))) 
		);

	return res;
}
// i j in "array notation" (change in row and collumn respectively)
double CalcTauYYCorrector(const std::vector<State>& field, size_t i, size_t j, size_t NPOINTSX, size_t NPOINTSY)
{
	const State& fieldij = field[j + (i * NPOINTSX)];
	const State& fieldijp = field[j + ((i-1) * NPOINTSX)];
	const State& fieldipj = field[j+1 + (i * NPOINTSX)];
	const State& fieldimj = field[j-1 + (i * NPOINTSX)];
	
	const double mu = CalcViscSutherland(muRef, tempRef, fieldij.T);

	const double res = (2.0/3.0) * mu * ( 
		((2.0/dy)*(fieldijp.v - fieldij.v)) - (((fieldipj.u - fieldimj.u)/(2.0*dx))) 
		);

	return res;
}
// i j in "array notation" (change in row and collumn respectively)
double CalcTauXYCorrectorE(const std::vector<State>& field, size_t i, size_t j, size_t NPOINTSX, size_t NPOINTSY)
{
	const State& fieldij = field[j + (i * NPOINTSX)];
	const State& fieldijp = field[j + ((i-1) * NPOINTSX)];
	const State& fieldipj = field[j+1 + (i * NPOINTSX)];
	const State& fieldijm = field[j + ((i+1) * NPOINTSX)];
	
	const double mu = CalcViscSutherland(muRef, tempRef, fieldij.T);

	const double res = mu * ( 
		((fieldipj.v - fieldij.v)/dx) + (((fieldijp.u - fieldijm.u)/(2.0*dy))) 
		);

	return res;
}
// i j in "array notation" (change in row and collumn respectively)
double CalcTauXYCorrectorF(const std::vector<State>& field, size_t i, size_t j, size_t NPOINTSX, size_t NPOINTSY)
{
	const State& fieldij = field[j + (i * NPOINTSX)];
	const State& fieldijp = field[j + ((i-1) * NPOINTSX)];
	const State& fieldipj = field[j+1 + (i * NPOINTSX)];
	const State& fieldimj = field[j-1 + (i * NPOINTSX)];
	
	const double mu = CalcViscSutherland(muRef, tempRef, fieldij.T);

	const double res = mu * ( 
		((fieldijp.u - fieldij.u)/dy) + (((fieldipj.v - fieldimj.v)/(2.0*dx))) 
		);

	return res;
}
// i j in "array notation" (change in row and collumn respectively)
double CalcQxCorrector(const std::vector<State>& field, size_t i, size_t j, size_t NPOINTSX, size_t NPOINTSY)
{
	const State& fieldij = field[j + (i * NPOINTSX)];
	const State& fieldipj = field[j+1 + (i * NPOINTSX)];
	
	const double mu = CalcViscSutherland(muRef, tempRef, fieldij.T);
	const double k = Calck(mu);
	
	const double res = -k * (fieldipj.T - fieldij.T)/dx;

	return res;
}
// i j in "array notation" (change in row and collumn respectively)
double CalcQyCorrector(const std::vector<State>& field, size_t i, size_t j, size_t NPOINTSX, size_t NPOINTSY)
{
	const State& fieldij = field[j + (i * NPOINTSX)];
	const State& fieldijp = field[j + ((i-1) * NPOINTSX)];
	
	const double mu = CalcViscSutherland(muRef, tempRef, fieldij.T);
	const double k = Calck(mu);
	
	const double res = -k * (fieldijp.T - fieldij.T)/dy;

	return res;
}

State CalcState(const State& U)
{
	State estado;

	estado.u = U.v/U.u;
	estado.v = U.P/U.u;

	const double e = (U.T / U.u) - (((estado.u*estado.u)+(estado.v+estado.v))/2.0);
	estado.T = e/Cv;
	estado.P = U.u * R * estado.T;

	return estado;
}

State CalcU(const State& estado)
{
	State U;
	const double Vquad = (estado.u*estado.u)+(estado.v+estado.v); 
	U.u = estado.P/(R*estado.T);
	U.v = U.u * estado.u;
	U.P = U.u* estado.v;
	U.T = U.u*((Cv * estado.T)+ (Vquad/2.0));

	return U;
}