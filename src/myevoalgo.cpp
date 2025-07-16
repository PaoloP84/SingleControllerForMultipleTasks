#include "myevoalgo.h"
#include "evonet.h"
#include <cmath>
#include <QTextStream>
#include <QFile>
#include <../unsupported/Eigen/MatrixFunctions>
#ifdef _WIN32
#include <Windows.h>
#else
#include <unistd.h>
#endif

using namespace std;

struct FitnessAndId
{
    /**
     * \brief Fitness
     */
    double fitness;

    /**
     * \brief Id
     */
    int id;
};

/**
 * \brief Lesser-than operator overloading for FitnessAndId
 */
bool operator<(FitnessAndId first, FitnessAndId second)
{
	return (first.fitness < second.fitness);
}

//------------------------------------- MYEVOALGO CLASS ------------------------------------------------------
//---------------------------------- Inherited from EVOGA ----------------------------------------------------

MyEvoAlgo::MyEvoAlgo()
	: farsa::Evoga()
	, m_exp(NULL)
	, m_nevalsteps(1000)
	, m_pop()
	, m_mutationRate(0.1)
	, m_crossoverRate(0.0)
	, m_range(1.0)
	, m_algo(0)
	, m_decay(false)
	, m_bestGenBuf()
	, m_statBuf()
	, m_evalstatBuf()
	, m_saveBestGenInfo(false)
	, m_bestGenInterval(1)
	, m_fitNoise(0.0)
	, m_refine(false)
	, m_numRefineIterations(0)
	, m_bests()
	, m_cartPole()
	, m_grid()
	, m_arena()
	, m_game()
	, m_rng(1)
	, m_rngFunct(1)
	, m_test(false)
{
}

MyEvoAlgo::~MyEvoAlgo()
{
}

void MyEvoAlgo::configure(farsa::ConfigurationParameters& params, QString prefix)
{
	// Calling parent function
	farsa::Evoga::configure(params, prefix);
    	m_nevalsteps = farsa::ConfigurationHelper::getInt(params, prefix + "nevalsteps", m_nevalsteps);
	m_mutationRate = farsa::ConfigurationHelper::getDouble(params, prefix + "mutationRate", m_mutationRate);
	if (m_mutationRate > 1.0)
	{
		m_mutationRate /= 100.0;
	}
	m_crossoverRate = farsa::ConfigurationHelper::getDouble(params, prefix + "crossoverRate", m_crossoverRate);
	if (m_crossoverRate > 1.0)
	{
		m_crossoverRate /= 100.0;
	}
	m_algo = farsa::ConfigurationHelper::getInt(params, prefix + "algo", m_algo);
	m_decay = farsa::ConfigurationHelper::getBool(params, prefix + "decay", m_decay);
	m_saveBestGenInfo = farsa::ConfigurationHelper::getBool(params, prefix + "saveBestGenInfo", m_saveBestGenInfo);
	if (m_saveBestGenInfo)
	{
		m_bestGenInterval = farsa::ConfigurationHelper::getInt(params, prefix + "bestGenInterval", m_bestGenInterval);
	}
	m_fitNoise = farsa::ConfigurationHelper::getDouble(params, prefix + "fitNoise", m_fitNoise);
	m_refine = farsa::ConfigurationHelper::getBool(params, prefix + "refine", m_refine);
	if (m_refine)
	{
		m_numRefineIterations = farsa::ConfigurationHelper::getInt(params, prefix + "numRefineIterations", m_numRefineIterations);
	}
	m_test = farsa::ConfigurationHelper::getBool(params, prefix + "test", m_test);
}

void MyEvoAlgo::describe(QString type)
{
	// Calling parent function
	farsa::Evoga::describe(type);

	Descriptor d = addTypeDescription(type, "A modified version of the steadyState algorithm that uses floating-point genes");
    	d.describeReal("mutationRate").def(0.1).help("The mutation rate");
	d.describeReal("crossoverRate").def(0.0).help("The crossover rate");
	d.describeInt("algo").limits(0,1).def(0).help("Algorithm id");
	d.describeBool("saveBestGenInfo").def(false).help("Flags whether or not the information about the best individual must be saved");
	d.describeInt("bestGenInterval").limits(1,INT_MAX).def(1).help("Interval for saving information about the current best individual");
	d.describeReal("fitNoise").def(0.0).help("The noise applied to the individual's fitness");
	d.describeBool("refine").def(false).help("Flags whether or not the refinement must be run");
	d.describeInt("numRefineIterations").limits(0,INT_MAX).def(0).help("Number of refinement iterations");
	d.describeBool("test").def(false).help("Test functions");
}

void MyEvoAlgo::postConfigureInitialization()
{
	// Calling parent function
	farsa::Evoga::postConfigureInitialization();
	// Cast from EvorobotExperiment to specific experiment class
	m_exp = dynamic_cast<farsa::CartPoleExperiment*>(exp);
	// Initialize cart pole
	m_cartPole.nmarkov_long = false;
	m_cartPole.initializeCartPole(false, true);
	// Resize population
	m_pop.resize(popSize * 2);
	for (int i = 0; i < (popSize * 2); i++)
		m_pop[i].resize(glen);
	m_bests.resize(nreproducing);
	for (int i = 0; i < nreproducing; i++)
		m_bests[i].resize(glen);

	farsa::ResourcesLocker locker(this);
	farsa::Evonet* evonet = getResource<farsa::Evonet>("evonet");
	m_range = evonet->getWrange();
}

void MyEvoAlgo::evolveAllReplicas()
{
	m_bestGenBuf.clear();
	m_statBuf.clear();
	m_fitBuf.clear();
	
	stopEvolution = false;
	// It doesn't matter the option chosen in evolutionType,
	// this class will always call my evolutionary Algorithm.
	if (m_test)
		test();
	else
	{
		if (m_algo == 0)
			evolveSSS();
		else if (m_algo == 1)
			evolveGGA();
		else if (m_algo == 2)
			evolveOpenAI();
		else
			evolveNES2();
	}
}

void MyEvoAlgo::setSeed(const int s)
{
	m_rng.setSeed(s);
	currentSeed = s;
}

farsa::RandomGenerator* MyEvoAlgo::getRng()
{
	return &m_rng;
}

farsa::RandomGenerator* MyEvoAlgo::getRngFunct()
{
	return &m_rngFunct;
}

void MyEvoAlgo::initPop()
{
	for (int i = 0; i < m_pop.size(); i++)
	{
		for (int j = 0; j < glen; j++)
		{
			m_pop[i][j] = getRng()->getDouble(-m_range,m_range);
		}
	}
}

void MyEvoAlgo::resetPop()
{
	for (int i = 0; i < m_pop.size(); i++)
	{
		m_pop[i].resize(0);
	}
	m_pop.resize(0);
}

void MyEvoAlgo::performMutation(int from, int to, int mut)
{
	int num;
	QVector<int> tmp(glen);
	// If <mut> is 0, no mutations are performed (i.e., simple copy)
    	if (mut == 0)
	{
		for (int i = 0; i < glen; i++) 
		{
			m_pop[to][i] = m_pop[from][i];
		}
    	}
    	else
	{
		if (getRng()->getDouble(0.0,1.0) < m_crossoverRate)
		{
			// Crossover
			num = getRng()->getInt(0, glen - 1);
			for (int i = num; i < glen; i++)
			{
				tmp[i - num] = m_pop[from][i];
			}
			for (int i = 0; i < num; i++)
			{
				tmp[glen - num + i] = m_pop[from][i];
			}
			// Copy with mutations
			for (int i = 0; i < glen; i++)
			{
				// Copy
				m_pop[to][i] = tmp[i];
				if (getRng()->getDouble(0.0,1.0) < m_mutationRate)
				{
					// Mutate
					m_pop[to][i] = getRng()->getDouble(-m_range,m_range);
				}
			}
		}
		else
		{
			// Mutation
			for (int i = 0; i < glen; i++)
			{
				if (getRng()->getDouble(0.0,1.0) < m_mutationRate)
				{
					m_pop[to][i] = getRng()->getDouble(-m_range,m_range);
		    		}
				else
				{
					m_pop[to][i] = m_pop[from][i];
				}
			}
		}
    	}
}

void MyEvoAlgo::immediateMutation(const QVector<float> from, QVector<float>& to, int mut)
{
	bool ok;
	int num;
	QVector<int> tmp(glen);
	// If <mut> is 0, no mutations are performed (i.e., simple copy)
    	if (mut == 0)
	{
		for (int i = 0; i < glen; i++) 
		{
			to[i] = from[i];
		}
    	}
    	else
	{
		if (getRng()->getDouble(0.0,1.0) < m_crossoverRate)
		{
			// Crossover
			num = getRng()->getInt(0, glen - 1);
			for (int i = num; i < glen; i++)
			{
				tmp[i - num] = from[i];
			}
			for (int i = 0; i < num; i++)
			{
				tmp[glen - num + i] = from[i];
			}
			// Copy with mutations
			for (int i = 0; i < glen; i++)
			{
				// Copy
				to[i] = tmp[i];
				if (getRng()->getDouble(0.0,1.0) < m_mutationRate)
				{
					// Mutate
					to[i] = getRng()->getDouble(-m_range,m_range);
				}
			}
		}
		else
		{
			// Mutation
			for (int i = 0; i < glen; i++)
			{
				if (getRng()->getDouble(0.0,1.0) < m_mutationRate)
				{
					to[i] = getRng()->getDouble(-m_range,m_range);
		    		}
				else
				{
					to[i] = from[i];
				}
			}
		}
    	}
}

void MyEvoAlgo::neutralMutation(const QVector<float> from, QVector<float>& to)
{
	// Extract a gene to be mutated (not an output gene)
	bool found = false;
	int geneIdx = -1;
	while (!found)
	{
		geneIdx = getRng()->getInt(0, glen);
		// Security check (perhaps it is useless!!!)
		found = (geneIdx < glen);
	}
	// Store the old value of the gene (to avoid the generation of the same individual)
	const int oldVal = from[geneIdx];
	// Modify the gene
	int currVal;
	bool modified = false;
	
	currVal = -1;
	while (!modified)
	{
		currVal = getRng()->getDouble(-m_range, m_range);
		modified = (currVal != oldVal);
	}
	for (int g = 0; g < glen; g++)
	{
		if (g == geneIdx)
		{
			to[g] = currVal;
		}
		else
		{
			to[g] = from[g];
		}
	}
}

double MyEvoAlgo::rastrigin(QVector<float> x)
{
	double f_out = 10.0 * glen;
	for (int i = 0; i < glen; i++)
		f_out += ((x[i] * x[i]) - (10.0 * cos(2.0 * PI_GRECO * x[i])));
	return f_out;
}

double MyEvoAlgo::rosenbrock(QVector<float> x)
{
	double f_out = 0.0;
	for (int i = 0; i < glen - 1; i++)
	{
		f_out += (100.0 * ((x[i + 1] - (x[i] * x[i])) * (x[i + 1] - (x[i] * x[i]))) + ((1 - x[i]) * (1 - x[i])));
	}
	return f_out;
}

double MyEvoAlgo::sphere(QVector<float> x)
{
	double f_out = 0.0;
	for (int i = 0; i < glen; i++)
		f_out += (x[i] * x[i]);
	return f_out;
}

double MyEvoAlgo::griewank(QVector<float> x)
{
	double f_out = 0.0;
	double f_sum = 0.0;
	double f_prod = 1.0;
	for (int i = 0; i < glen; i++)
	{
		f_sum += (x[i] * x[i]);
		f_prod *= cos(x[i] / sqrt(i + 1));
	}
	f_out = 1.0 + f_sum / 4000.0 - f_prod;
	return f_out;
}

double MyEvoAlgo::ackley(QVector<float> x)
{
	const double a = 20.0;
	const double b = 0.2;
	const double c = 2.0 * PI_GRECO;
	double f_out = 0.0;
	double f_sum = 0.0;
	double f_sum_2 = 0.0;
	for (int i = 0; i < glen; i++)
	{
		f_sum += (x[i] * x[i]);
		f_sum_2 += cos(c * x[i]);
	}
	f_out = -a * std::exp(-b * sqrt(f_sum / glen)) - std::exp(f_sum_2 / glen) + a + std::exp(1);
	return f_out;
}

double MyEvoAlgo::evaluateParity(const QVector<float> ind, const QVector< QVector<int> > bitStrings, int& steps)
{
    	double fitness = 0.0;
	const int numInputs = bitStrings.size();
	int inputSize;
	double output;
	double expOut;
	int count;
	farsa::ResourcesLocker locker(this);
	farsa::Evonet* evonet = getResource<farsa::Evonet>("evonet");
	// Set parameters to network
	evonet->setParameters(ind);
	steps = 0;
	for (int i = 0; i < numInputs; i++)
	{
		evonet->resetNet();
		inputSize = bitStrings[i].size();
		// Set the input
		for (int j = 0; j < inputSize; j++)
			evonet->setInput(j, bitStrings[i][j]);
		// Activate the net
		evonet->updateNet();
		// Get the output
		output = evonet->getOutput(0);
		if (output >= 0.0)
			output = 1.0;
		else
			output = 0.0;
		// Compute the expected out for the current input
		for (int j = 0; j < inputSize; j++)
		{
			if (bitStrings[i][j] == 1)
				count++;
		}
		expOut = ((count % 2) == 0) ? 1 : 0;
		fitness += (expOut == output) ? 1.0 : 0.0;
		steps++;
	}
	//fitness /= numInputs;
    	return numInputs - fitness;//(1.0 - fitness); // Minimization
}

double MyEvoAlgo::evaluateDpole(const QVector<float> ind, int& steps)
{
	double fitness = 0.0;
	int i;
	const double len = m_cartPole.getTaskLength();
	double fit;
	farsa::ResourcesLocker locker(this);
	farsa::Evonet* evonet = getResource<farsa::Evonet>("evonet");
	// Set parameters to network
	evonet->setParameters(ind);
	steps = 0;
	// Evaluate network
	for (i = 0; i < DPOLE_EPS; i++)
	{
		evonet->resetNet();
		fit = m_cartPole.evalNet(evonet, i);
		fitness += (len - fit);//1.0 - (fit / len); // Minimization [0,1000]
		steps += int(fit);
	}
	fitness /= DPOLE_EPS;
    	return fitness;
}

double MyEvoAlgo::evaluateGrid(const QVector<float> ind, int& steps)
{
	double fitness = 0.0;
	int i;
	double fit;
	double dist;
	farsa::ResourcesLocker locker(this);
	farsa::Evonet* evonet = getResource<farsa::Evonet>("evonet");
	// Set parameters to network
	evonet->setParameters(ind);
	//m_randGrid.setSeed(LOCAL_SEED);
	steps = 0;
	// Evaluate network
	for (i = 0; i < GRID_EPS; i++)
	{
		evonet->resetNet();
		fit = m_grid.evalNet(evonet, i, dist);
		fitness += dist;
		steps += int(fit);
	}
	fitness /= GRID_EPS;
    	return fitness; // Minimization of the distance
}

double MyEvoAlgo::evaluateFunct(const QVector<float> ind, int& steps)
{
	const int size = ind.size();
	double fitness = 0.0;
	double fit_ra, fit_ro, fit_sp, fit_gr, fit_ac;
	steps = 0;
	// Rastrigin
	fit_ra = rastrigin(ind);
	fitness += fit_ra;
	steps++;
	// Rosenbrock
	fit_ro = rosenbrock(ind);
	fitness += fit_ro;
	steps++;
	// Sphere
	fit_sp = sphere(ind);
	fitness += fit_sp;
	steps++;
	// Griewank
	fit_gr = griewank(ind);
    	fitness += fit_gr;
	steps++;
	// Ackley
	fit_ac = ackley(ind);
    	fitness += fit_ac;
	steps++;
	//printf("%lf\t%lf\t%lf\t%lf\t%lf\n", fit_ra, fit_ro, fit_sp, fit_gr, fit_ac);
	fitness /= 5;
	return fitness;
}

double MyEvoAlgo::evaluateExploration(const QVector<float> ind, int& steps)
{
	/*double fitness = 0.0;
	int i;
	double fit;
	int num;
	farsa::ResourcesLocker locker(this);
	farsa::Evonet* evonet = getResource<farsa::Evonet>("evonet");
	// Set parameters to network
	evonet->setParameters(ind);
	//m_randGrid.setSeed(LOCAL_SEED);
	steps = 0;
	// Evaluate network
	evonet->resetNet();
	fit = m_arena.evalNet(evonet, num);
	fitness = (double)num / (ARENA_SIZE * ARENA_SIZE);
	steps += int(fit);
    	return fitness;*/
}

double MyEvoAlgo::evaluateGame(const QVector<float> ind, int& steps)
{
	/*double fitness = 0.0;
	int i;
	double fit;
	farsa::ResourcesLocker locker(this);
	farsa::Evonet* evonet = getResource<farsa::Evonet>("evonet");
	// Set parameters to network
	evonet->setParameters(ind);
	//m_game.setSeed(LOCAL_SEED);
	steps = 0;
	// Evaluate network
	for (i = 0; i < GAME_EPS; i++)
	{
		evonet->resetNet();
		fitness += m_game.evalNet(evonet, i);
		steps += GAME_STEPS;
	}
	fitness /= GAME_EPS;
    	return fitness;*/
}

double MyEvoAlgo::evaluate(const QVector<float> ind, const QVector< QVector<int> > bitStrings, int& steps, QVector<float>& taskFit)
{
	double fitness = 0.0;
	int numSteps;

	if (taskFit.size() != NUM_TASKS)
		taskFit.resize(NUM_TASKS);

	// Compute tasks' fitnesses and update steps
	// Bit-parity
	double fit_parity = evaluateParity(ind, bitStrings, numSteps);
	steps += numSteps;
	// Double-pole balancing
	double fit_dpole = evaluateDpole(ind, numSteps);
	steps += numSteps;
	// Grid navigation
	double fit_grid = evaluateGrid(ind, numSteps);
	steps += numSteps;
	// Prediction
	double fit_funct = evaluateFunct(ind, numSteps);
	//printf("fit_funct = %d\n", fit_funct);
	steps += numSteps;
	//double fit_expl = evaluateExploration(ind, numSteps);
	//steps += numSteps;
	// Game play
	//double fit_game = evaluateGame(ind, numSteps);
	//steps += numSteps;
	// Compute overall fitness as the sum of the fitnesses
	fitness = (fit_parity + fit_dpole + fit_grid + fit_funct);// / NUM_TASKS;
	// Save detailed fitnesses
	taskFit[0] = fit_parity;
	taskFit[1] = fit_dpole;
	taskFit[2] = fit_grid;
	taskFit[3] = fit_funct;
	//taskFit[4] = fit_game;
	/*for (int i = 0; i < NUM_TASKS; i++)
		printf("fit[%d] = %lf\n", i, taskFit[i]);*/
	return fitness;
}

double MyEvoAlgo::refine(const QVector<float> original, double fitness, QVector<float> taskFit, const QVector< QVector<int> > bitStrings, QVector<float>& novel, int& steps, QVector<float>& novelTaskFit, bool& solved)
{
	double novelFit;
	double currFit;
	// Compute all the bit strings of length <m_input_size>
    	QVector<float> orig = original;
	QVector<float> evalInd(glen);
	QVector<float> evalTaskFit(NUM_TASKS);
	int iter = 0;
	double fit = 0.0;
	double storedFit = fitness;
	int nsteps;
	// Initialize task fit
	for (int i = 0; i < NUM_TASKS; i++)
		novelTaskFit[i] = taskFit[i];
	// Reset temporary individual
	invalidFill(evalInd);
	while (iter < m_numRefineIterations)
	{
		// We must generate new noise at each iteration to the original fitness
		//currFit = storedFit * (1.0 + getRng()->getDouble(-m_fitNoise, m_fitNoise));
		// Mutate the individual
		//immediateMutation(orig, evalInd, 1);
		neutralMutation(orig, evalInd);
		// Evaluate
		nsteps = 0;
		fit = evaluate(evalInd, bitStrings, nsteps, evalTaskFit);
		steps += nsteps;
		if ((fit == (double)NUM_TASKS) && !solved)
		{
			solved = true;
			// Save the number of evaluation steps
			char evalfname[1024];
			sprintf(evalfname, "evalS%d.txt", currentSeed);
			FILE* evalfp = fopen(evalfname, "w");
			if (evalfp != NULL)
			{
				fprintf(evalfp, "%d", steps);
				fclose(evalfp);
			}
		}
		// Add noise
		/*double noisyFit = fit * (1.0 + getRng()->getDouble(-m_fitNoise, m_fitNoise));
		if (noisyFit < 0.0)
			noisyFit = 0.0;
		if (noisyFit > 1.0)
			noisyFit = 1.0;
		if (noisyFit <= currFit)*/
		if (fit <= storedFit)
		{
			immediateMutation(evalInd, orig, 0);
			storedFit = fit;
			//currFit = noisyFit;
			for (int i = 0; i < NUM_TASKS; i++)
				novelTaskFit[i] = evalTaskFit[i];
		}
		iter++;
	}
	novel = orig;
	novelFit = storedFit;
	//printf("fit: %lf vs %lf\n", fitness, novelFit);
	return novelFit;
}

void MyEvoAlgo::loadGen(FILE *fp, int ind)
{
	float v;

	fscanf(fp, "DYNAMICAL NN\n");
	for (int g = 0; g < glen; g++) 
	{
		fscanf(fp, "%f\n", &v);
		m_pop[ind][g] = v;
	}
	fscanf(fp, "END\n");
}

int MyEvoAlgo::loadAllGen(int gen, char* filew)
{
	char filename[512];
	char message[512];
	char flag[512];
	FILE* fp;
	if (gen >= 0)
	{
		sprintf(filename, "G%dS%d.gen", gen, currentSeed);
	}
	else
	{
		sprintf(filename, "%s", filew);
	}

	fp = fopen(filename, "r");
	if (fp != NULL) 
	{
		resetPop();
		bool cond = true;
		while (cond) {
			flag[0] = '\0';
			fscanf(fp, "%s : %s\n", flag, message);
			if (strcmp(flag, "**NET") == 0) 
			{
				QVector<float> ind(glen);
				m_pop.append(ind);
				loadGen(fp, m_pop.size() - 1);
			} 
			else 
			{
				cond = false;
			}
		}
		farsa::Logger::info(QString("Loaded ind: %1").arg(m_pop.size()));
		fclose(fp);
	} 
	else 
	{
		farsa::Logger::error(QString("File %1 could not be opened").arg(filename));
	}

	loadedIndividuals = m_pop.size();

	return m_pop.size();
}

void MyEvoAlgo::loadBestInd(QVector<float>& ind)
{
	FILE* fp;
	char filename[64];
	int i;
	float v;

	if (ind.size() != glen)
		ind.resize(glen);

	sprintf(filename, "bestIndS%d.gen", currentSeed);
	fp = fopen(filename, "r");
	if (fp != NULL)
	{
		for (int i = 0; i < glen; i++)
		{
			fscanf(fp, "%f\n", &v);
			ind[i] = v;
		}
		fclose(fp);
	}
}

void MyEvoAlgo::saveBestInd(QVector<float> ind)
{
	FILE* fp;
	char filename[64];
	int i;
	sprintf(filename, "bestIndS%d.gen", currentSeed);
	fp = fopen(filename, "w");
	if (fp != NULL)
	{
		for (int i = 0; i < glen; i++)
			fprintf(fp, "%f\n", ind[i]);
		fclose(fp);
	}
}

void MyEvoAlgo::saveGen(FILE *fp, int ind)
{
    int j;
    fprintf(fp, "DYNAMICAL NN\n");
    for (j = 0; j < glen; j++)
	{
        fprintf(fp, "%f\n", m_pop[ind][j]);
	}
    fprintf(fp, "END\n");
}

void MyEvoAlgo::saveAllGen()
{
    FILE *fp;
    char filename[64];
    int i;

    sprintf(filename, "G%dS%d.gen", cgen, currentSeed);
    if ((fp = fopen(filename, "w+")) == NULL) 
	{
        farsa::Logger::error(QString("Cannot open file %1").arg(filename));
    } 
	else {
        //we save
        for(i = 0; i < popSize; i++) 
		{
            fprintf(fp, "**NET : %d_%d_%d.wts\n", cgen, 0, i);
            saveGen(fp, i);
        }
        fclose( fp );
    }
}

void MyEvoAlgo::computeFStat2()
{
	int i;
	double min, max, av;

	//min=max=tfitness[0]/ntfitness[0];
	//try to fix a problem
	min=999999999.00;
	max=-999999999.00;
	av=0.0;

	for(i=0;i<this->popSize;i++) {
		if((tfitness[i]/ntfitness[i])<min) {
			min=tfitness[i]/ntfitness[i];
		}
		if((tfitness[i]/ntfitness[i])>max) max=tfitness[i]/ntfitness[i];
		av+=(tfitness[i]/ntfitness[i]);
	}
	this->faverage=av/(double)this->popSize;
	this->fmax=max;
	this->fmin=min;
	this->statfit[this->cgen][0]=this->fmin;
	this->statfit[this->cgen][1]=this->faverage;
	this->statfit[this->cgen][2]=this->fmax;

	if (this->fmin < this->fbest) {
		this->fbest = this->fmin;
		this->fbestgen = this->cgen;
	}
}

void MyEvoAlgo::saveFStat()
{
	m_statBuf += QString("%1 %2 %3\n").arg(fmax).arg(faverage).arg(fmin);
}

void MyEvoAlgo::saveFitStats(QVector<float> fit)
{
	m_fitBuf += QString("%1 %2 %3 %4\n").arg(fit[0]).arg(fit[1]).arg(fit[2]).arg(fit[3]);
}

void MyEvoAlgo::saveEvalStats(int steps, double fit)
{
	m_evalstatBuf += QString("%1 %2\n").arg(steps).arg(fit);
}

void MyEvoAlgo::invalidFill(QVector<float>& ind)
{
	for (int i = 0; i < ind.size(); i++)
	{
		ind[i] = -1.0;
	}
}

/*
 * Main function of the Genetic Algorithm (Steady State Version with the possibility
 * of performing the Annealing routine...)
 */
void MyEvoAlgo::evolveSSS()
{
	int rp; //replication
	int gn; //generation
	int cstep; // evaluation step
    	int id; //individuals
	double fit;
	int startGeneration = 0;
	char statfile[128];
	char genFile[128];
	char filename[64];
    	int steps;
	bool solved;
	double bestfit;
	int bestid;
	QVector<float> bestTaskFit(NUM_TASKS);
	QVector<float> bestInd(glen);

	const int numInputs = (int)pow(2.0,(double)NUM_INPUTS);
	QVector<QVector<int> > bitStrings(numInputs);
	int** tmpStrings = generateBitstrings(NUM_INPUTS);
	for (int i = 0; i < numInputs; i++)
	{
		bitStrings[i].resize(NUM_INPUTS);
		for (int j = 0; j < NUM_INPUTS; j++)
		{
			bitStrings[i][j] = tmpStrings[i][j];
		}
	}

	QVector< QVector<float> > taskFit(popSize * 2);
	for (id = 0; id < popSize * 2; id++)
		taskFit[id].resize(NUM_TASKS);
	
	farsa::Logger::info("EVOLUTION: steady state with custom evaluation function");
    	farsa::Logger::info("Number of replications: " + QString::number(nreplications));

	// Individual to be evaluated
	QVector<int> evalInd(glen); // Vector to be used to avoid overwritings!!
	// replications
    	for(rp = 0; rp < nreplications; rp++) 
	{
		startGeneration = 0;
        	setSeed(getStartingSeed() + rp);
		farsa::Logger::info(QString("Replication %1, seed: %2").arg(rp + 1).arg(getStartingSeed() + rp));
		resetGenerationCounter();
		// Initialize the population (i.e. only the parents, offspring are mutated copies)
        	initPop();
		// Set fbest to a very low value
		this->fbest = 99999.0;

		emit startingReplication( rp );
       
		QTime evotimer;
        	
		// Reset tfitness and ntfitness
		for (int i = 0; i < m_pop.size(); i++) 
		{
		    tfitness[i] = 0.0;
		    ntfitness[i] = 0.0;
		}
		//code to recovery a previous evolution: Experimental
		sprintf(statfile, "statS%d.fit", getStartingSeed() + rp);
		//now check if the file exists
		farsa::DataChunk statTest(QString("stattest"), Qt::blue, 2000, false);
		if (statTest.loadRawData(QString(statfile),0))
		{
		    startGeneration = statTest.getIndex();
		    sprintf(genFile,"G%dS%d.gen", startGeneration, getStartingSeed() + rp);
		    farsa::Logger::info("Recovering from startGeneration: " + QString::number(startGeneration));
		    farsa::Logger::info(QString("Loading file: ") + genFile);
		    loadAllGen(-1, genFile);
		    cgen = startGeneration;
		    emit recoveredInterruptedEvolution( QString(statfile) );
		} //end evolution recovery code

		// Buffers to save statistics, genotypes, fitnesses and the like
		m_bestGenBuf = "";
		m_statBuf = "";
		m_fitBuf = "";
		m_evalstatBuf = "";
		QString tmpStatBuf = "";
		QString tmpFitBuf = "";
		QString tmpEvalStatBuf = "";

		// Set the cart pole
		m_cartPole.setCartPole(getStartingSeed() + rp);
		// generations
		gn = startGeneration;
		cstep = 0;
		bestfit = 999999999.0;
		bestid = -1;
		solved = false;
		while (cstep < m_nevalsteps)
        	//for(gn=startGeneration;gn<nogenerations;gn++) 
		{
			evotimer.restart();
		    	vector<int> identity;
		    	farsa::Logger::info(" Generation " + QString::number(gn + 1));
		    	exp->initGeneration(gn);
		    	if ( commitStep() ) 
			{ 
				return;
			}
			//individuals
			for(id = 0; id < popSize; id++) 
			{
				identity.push_back(id);
				performMutation(id, popSize + id, 1);
				// Test the individual and its child
				// Parents are tested only at the first generation
				if ((gn == startGeneration) || (m_fitNoise > 0.0))
				{
					// evaluate the parent
					steps = 0;
					fit = evaluate(m_pop[id], bitStrings, steps, taskFit[id]);
					// update its fitness
					tfitness[id] = fit;
					ntfitness[id] = 1.0;
					if (isStopped()) 
					{ // stop evolution
						return;
					}
					cstep += steps;
					if ((fit == (double)NUM_TASKS) && !solved)
					{
						solved = true;
						// Save the number of evaluation steps
						char evalfname[1024];
						sprintf(evalfname, "evalS%d.txt", currentSeed);
						FILE* evalfp = fopen(evalfname, "w");
						if (evalfp != NULL)
						{
							fprintf(evalfp, "%d", cstep);
							fclose(evalfp);
						}
					}
					// Check whether or not individual is better than current best
					if (fit < bestfit)
					{
						bestfit = fit;
						bestid = id;
						for (int i = 0; i < NUM_TASKS; i++)
							bestTaskFit[i] = taskFit[id][i];
						for (int i = 0; i < glen; i++)
							bestInd[i] = m_pop[id][i];
					}
				}
				// evaluate the child
				steps = 0;
				fit = evaluate(m_pop[popSize + id], bitStrings, steps, taskFit[popSize + id]);
				// update its fitness
				tfitness[popSize + id] = fit;
				ntfitness[popSize + id] = 1.0;
				if (isStopped()) 
				{ // stop evolution
					return;
				}
				cstep += steps;
				if ((fit == (double)NUM_TASKS) && !solved)
				{
					solved = true;
					// Save the number of evaluation steps
					char evalfname[1024];
					sprintf(evalfname, "evalS%d.txt", currentSeed);
					FILE* evalfp = fopen(evalfname, "w");
					if (evalfp != NULL)
					{
						fprintf(evalfp, "%d", cstep);
						fclose(evalfp);
					}
				}
				// Check whether or not individual is better than current best
				if (fit < bestfit)
				{
					bestfit = fit;
					bestid = (popSize + id);
					for (int i = 0; i < NUM_TASKS; i++)
						bestTaskFit[i] = taskFit[popSize + id][i];
					for (int i = 0; i < glen; i++)
						bestInd[i] = m_pop[popSize + id][i];
				}
			}
            		exp->endGeneration(gn);
            		if ( commitStep() ) 
			{ 
				return; 
			}
			

			// ========= Selection part ========== //
		    	// Finally, we look for the worst individuals (parents) and substitute them with best children.
		    	// What we do is: we order both the parents and the children in descending order, then we take the
		    	// popSize best individuals. This is not the same as what is done in the sequential version of
		    	// the algorithm, but should be similar. We overwrite the worst parents with the best children
			random_shuffle(identity.begin(),identity.end());
				
		    	QVector<FitnessAndId> parents(popSize);
		    	QVector<FitnessAndId> children(popSize);
		    	for (int i = 0; i < popSize; i++)
			{
				parents[i].fitness = tfitness[identity[i]] / ntfitness[identity[i]];
                		parents[i].id = identity[i];
                		parents[i].fitness *= (1.0 + getRng()->getDouble(-m_fitNoise, m_fitNoise));
            		}

            		for (int i = 0; i < popSize; i++)
			{
				children[i].fitness = tfitness[popSize + identity[i]] / ntfitness[popSize + identity[i]];
                		children[i].id = popSize + identity[i];
                		children[i].fitness *= (1.0 + getRng()->getDouble(-m_fitNoise, m_fitNoise));
            		}
			// Sorting both parents and children. They are sorted in ascending order but we need the best
            		// individuals (those with the highest fitness) first
            		qSort(parents);   //Error de operadores aqui#include <QTime>
            		qSort(children);
			// Store the list of selected and replaced individuals
			QVector<int> selected(popSize);
			QVector<float> selectedFit(popSize);
			QVector<int> replaced(popSize);
            		int p = 0;//popSize - 1;
            		int c = 0;//popSize - 1;
            		for (int i = 0; i < popSize; i++)
			{
                		if (parents[p].fitness < children[c].fitness) 
				{
                    			selected[i] = parents[p].id;
					selectedFit[i] = parents[p].fitness;
					replaced[i] = parents[p].id;
					p++;
            			} 
				else 
				{
					selected[i] = children[c].id;
					selectedFit[i] = children[c].fitness;
					replaced[i] = parents[popSize - 1 - c].id;
                    			c++;
                		}
            		}
			// Check whether or not to perform refinement
			if (m_refine)
			{
				// Refinement
				QVector< QVector<float> > refineTaskFit(popSize);
				for (id = 0; id < popSize; id++)
					refineTaskFit[id].resize(NUM_TASKS);
				// Parents
				for (id = 0; id < popSize; id++)
				{
					// Get the individual index
					const int index = selected[id];
					QVector<float> original(glen);
					immediateMutation(m_pop[index], original, 0);
					double fitness = tfitness[index];
					QVector<float> novel(glen);
					double novelFit = refine(original, fitness, taskFit[index], bitStrings, novel, cstep, refineTaskFit[id], solved);
					// Update the individual
					immediateMutation(novel, m_pop[index], 0);
					// Update structure data
					selectedFit[id] = novelFit;
					tfitness[index] = selectedFit[id];
					for (int i = 0; i < NUM_TASKS; i++)
						taskFit[index][i] = refineTaskFit[id][i];
					// Check whether or not the individual returned by the annealing routine is better than current best
					if (novelFit < bestfit)
					{
						bestfit = novelFit;
						bestid = index;
						for (int i = 0; i < NUM_TASKS; i++)
							bestTaskFit[i] = taskFit[index][i];//refineTaskFit[id][i];
						for (int i = 0; i < glen; i++)
							bestInd[i] = m_pop[index][i];
					}
				}
			}
			// Swap individuals
            		for (id = 0; id < popSize; id++)
			{
                		performMutation(selected[id], replaced[id], 0);
                		tfitness[replaced[id]] = tfitness[selected[id]];
                		ntfitness[replaced[id]] = ntfitness[selected[id]];
				for (int i = 0; i < NUM_TASKS; i++)
					taskFit[replaced[id]][i] = taskFit[selected[id]][i];
            		}
			// Compute and save statistics
		    	computeFStat2();
		    	saveFStat();

            		emit endGeneration( cgen, fmax, faverage, fmin );
            		if (commitStep()) 
			{
                		return; // stop the evolution process
            		}

            		cgen++;

            		farsa::Logger::info(QString("Generation %1 took %2 minutes - Generation's best fitness %3 - Best fitness = %4 - Evaluation steps = %5").arg(gn+1).arg((double)evotimer.elapsed()/60000.0, 0, 'f', 2).arg(fmin).arg(bestfit).arg(cstep));
            		fflush(stdout);

			saveEvalStats(cstep, bestfit);

			// Save detailed stats about task fitnesses
			saveFitStats(bestTaskFit);
 

			// Check whether the task is finished or not
			if (cstep >= m_nevalsteps)
			{
				if (!solved)
				{
					// Save the number of evaluation steps
					char evalfname[1024];
					sprintf(evalfname, "evalS%d.txt", currentSeed);
					FILE* evalfp = fopen(evalfname, "w");
					if (evalfp != NULL)
					{
						fprintf(evalfp, "%d", cstep);
						fclose(evalfp);
					}
					else
					{
						farsa::Logger::error("ERROR IN OPENING FILE " + QString(evalfname));
					}
				}
				
				tmpStatBuf += m_statBuf;
				// Saving files
				{
					const QString statFilename = QString("statS%1.fit").arg(currentSeed);
					QFile statFile(statFilename);
					if (!statFile.open(QIODevice::WriteOnly | QIODevice::Append | QIODevice::Truncate)) {
						farsa::Logger::warning("Cannot save stats into " + statFilename);
					} else {
						QTextStream s(&statFile);
						s << tmpStatBuf;
					}
					m_statBuf.clear();
				}

				{
					tmpFitBuf += m_fitBuf;
					// And best generalization fitness
					const QString fitFilename = QString("taskFitS%1.fit").arg(currentSeed);
					QFile fitFile(fitFilename);
					if (!fitFile.open(QIODevice::WriteOnly | QIODevice::Append | QIODevice::Truncate)) {
						farsa::Logger::warning("Cannot save stats into " + fitFilename);
					} else {
						QTextStream s(&fitFile);
						s << tmpFitBuf;
					}
					m_fitBuf.clear();
				}

				tmpEvalStatBuf += m_evalstatBuf;
				// Saving files
				{
					const QString evalStatFilename = QString("fitS%1.txt").arg(currentSeed);
					QFile evalStatFile(evalStatFilename);
					if (!evalStatFile.open(QIODevice::WriteOnly | QIODevice::Append | QIODevice::Truncate)) {
						farsa::Logger::warning("Cannot save stats into " + evalStatFilename);
					} else {
						QTextStream s(&evalStatFile);
						s << tmpEvalStatBuf;
					}
					//if (m_resetBufNumGen > 1)
					{
						m_evalstatBuf.clear();
					}
				}
			}
			
			// Saving files
			if (m_saveBestGenInfo)
			{
				if ((gn == 0) || (((gn + 1) % m_bestGenInterval) == 0))
				{
					const QString bestGenFilename = QString("B0S%1G%2.gen").arg(currentSeed).arg(gn + 1);
					QFile bestGenFile(bestGenFilename);
					if (!bestGenFile.open(QIODevice::WriteOnly | QIODevice::Append | QIODevice::Truncate)) {
						farsa::Logger::warning("Cannot save best genomes into " + bestGenFilename);
					} else {
						QTextStream s(&bestGenFile);
						s << m_bestGenBuf;
					}
					m_bestGenBuf.clear();
				}
			}
			gn++;
        	}
		saveBestInd(bestInd);
        	saveAllGen();
    	}
}

//Reproduce individuals with higher ranking
void MyEvoAlgo::reproduce()
{
	//to do
	//selecting best fathers
	int i;
	int bi,bx;
	double bn;
	int idx;
	int g;
	double coin;
	int point;

	char sbuffer[64];
	FILE *fp;

	// Vector storing indices
	vector<int> identity;

	//first of all we compute fitness stat
	this->computeFStat2();

	// Fill identity vector
	for(i = 0; i < this->popSize; i++) 
		identity.push_back(i); // Indices are stored in ascending order
	// Randomize indices to avoid selecting always the same agents based on their indices
	random_shuffle(identity.begin(),identity.end()); // This operation changes the order of indices in the array

	// Select best <nreproducing> individuals as fathers
	for(bi = 0; bi < this->nreproducing; bi++)
	{
		bn = 999999999.0;
		bx = -1; // individual to be copied
		for(i = 0; i < this->popSize; i++) 
		{
			idx = identity[i];
			if(tfitness[idx] / ntfitness[idx] < bn) 
			{
				bn = tfitness[idx] / ntfitness[idx];
				bx = idx;
			}
		}

		for (g = 0; g < glen; g++)
			m_bests[bi][g] = m_pop[bx][g];

		tfitness[bx] = 999999999.0;
	}

	//reproducing best
	bx = 0;
	for(bi = 0; bi < this->nreproducing; bi++)
	{
		for(i = 0; i < this->noffspring; i++) 
		{
			if (i == 0)
			{
				// Elitism
				for (g = 0; g < glen; g++)
					m_pop[bx][g] = m_bests[bi][g];
			}
			else
			{
				coin = getRng()->getDouble(0.0,1.0);
				if (coin < m_crossoverRate)
				{
					// Apply crossover
					point = getRng()->getInt(0, glen - 1);
					for (g = point; g < glen; g++)
						m_pop[bx][g - point] = m_bests[bi][g];
					for (g = 0; g < point; g++)
						m_pop[bx][glen - point + g] = m_bests[bi][g];
					// Try to apply mutations too
					for (g = 0; g < glen; g++)
					{
						coin = getRng()->getDouble(0.0,1.0);
						if (coin < m_mutationRate)
							// Apply mutation
							m_pop[bx][g] = getRng()->getDouble(-m_range, m_range);
					}
				}
				else
				{
					for (g = 0; g < glen; g++)
					{
						m_pop[bx][g] = m_bests[bi][g];
						coin = getRng()->getDouble(0.0,1.0);
						if (coin < m_mutationRate)
							// Apply mutation
							m_pop[bx][g] = getRng()->getDouble(-m_range, m_range);
					}
				}
				//this->getGenome(bi, bx, 1);
			}
			bx++;
		}
	}
	// resetting fitness
	for (i = 0; i < this->popSize; i++)
	{
		tfitness[i] = 0.0;
		ntfitness[i] = 0.0;
	}
	// Save stats
	this->saveFStat();
	// Update generation counter
	cgen++;
}

/*
 * Main function of the Genetic Algorithm (Steady State Version with the possibility
 * of performing the Annealing routine...)
 */
void MyEvoAlgo::evolveGGA()
{
	int rp; //replication
	int gn; //generation
	int cstep; // evaluation step
    	int id; //individuals
	double fit;
	int startGeneration = 0;
	char statfile[128];
	char genFile[128];
	char filename[64];
	FILE* fp;
    	int steps;
	bool solved;
	double bestfit;
	int bestid;
	QVector<float> bestTaskFit(NUM_TASKS);
	QVector<float> bestInd(glen);

	// Resizing population
	m_pop.resize(popSize);

	const int numInputs = (int)pow(2.0,(double)NUM_INPUTS);
	QVector<QVector<int> > bitStrings(numInputs);
	int** tmpStrings = generateBitstrings(NUM_INPUTS);
	for (int i = 0; i < numInputs; i++)
	{
		bitStrings[i].resize(NUM_INPUTS);
		for (int j = 0; j < NUM_INPUTS; j++)
		{
			bitStrings[i][j] = tmpStrings[i][j];
		}
	}

	QVector< QVector<float> > taskFit(popSize);
	for (id = 0; id < popSize; id++)
		taskFit[id].resize(NUM_TASKS);
	
	farsa::Logger::info("EVOLUTION: GGA with custom evaluation function");
    	farsa::Logger::info("Number of replications: " + QString::number(nreplications));

	// Individual to be evaluated
	QVector<int> evalInd(glen); // Vector to be used to avoid overwritings!!
	// replications
    	for(rp = 0; rp < nreplications; rp++) 
	{
		startGeneration = 0;
        	setSeed(getStartingSeed() + rp);
		farsa::Logger::info(QString("Replication %1, seed: %2").arg(rp + 1).arg(getStartingSeed() + rp));
		resetGenerationCounter();
		// Initialize the population (i.e. only the parents, offspring are mutated copies)
        	initPop();
		// Set fbest to a very low value
		this->fbest = 99999999999.0;

		emit startingReplication( rp );
       
		QTime evotimer;
        	
		// Reset tfitness and ntfitness
		for (int i = 0; i < m_pop.size(); i++) 
		{
		    tfitness[i] = 0.0;
		    ntfitness[i] = 0.0;
		}
		//code to recovery a previous evolution: Experimental
		sprintf(statfile, "statS%d.fit", getStartingSeed() + rp);
		//now check if the file exists
		farsa::DataChunk statTest(QString("stattest"), Qt::blue, 2000, false);
		if (statTest.loadRawData(QString(statfile),0))
		{
		    startGeneration = statTest.getIndex();
		    sprintf(genFile,"G%dS%d.gen", startGeneration, getStartingSeed() + rp);
		    farsa::Logger::info("Recovering from startGeneration: " + QString::number(startGeneration));
		    farsa::Logger::info(QString("Loading file: ") + genFile);
		    loadAllGen(-1, genFile);
		    cgen = startGeneration;
		    emit recoveredInterruptedEvolution( QString(statfile) );
		} //end evolution recovery code

		// Buffers to save statistics, genotypes, fitnesses and the like
		m_bestGenBuf = "";
		m_statBuf = "";
		m_fitBuf = "";
		m_evalstatBuf = "";
		QString tmpStatBuf = "";
		QString tmpFitBuf = "";
		QString tmpEvalStatBuf = "";

		// Set the cart pole
		m_cartPole.setCartPole(getStartingSeed() + rp);
		// generations
		gn = startGeneration;
		cstep = 0;
		bestfit = 999999999.0;
		bestid = -1;
		solved = false;
		while (cstep < m_nevalsteps)
        	//for(gn=startGeneration;gn<nogenerations;gn++) 
		{
			evotimer.restart();
		    	vector<int> identity;
		    	farsa::Logger::info(" Generation " + QString::number(gn + 1));
		    	exp->initGeneration(gn);
		    	if ( commitStep() ) 
			{ 
				return;
			}
			//individuals
			for(id = 0; id < popSize; id++) 
			{
				// evaluate the parent
				steps = 0;
				fit = evaluate(m_pop[id], bitStrings, steps, taskFit[id]);
				//printf("gen %d id %d fit %lf\n", cgen, id, fit);
				// update its fitness
				tfitness[id] = fit;
				ntfitness[id] = 1.0;
				if (isStopped()) 
				{ // stop evolution
					return;
				}
				cstep += steps;
				if ((fit == (double)NUM_TASKS) && !solved)
				{
					solved = true;
					// Save the number of evaluation steps
					char evalfname[1024];
					sprintf(evalfname, "evalS%d.txt", currentSeed);
					FILE* evalfp = fopen(evalfname, "w");
					if (evalfp != NULL)
					{
						fprintf(evalfp, "%d", cstep);
						fclose(evalfp);
					}
				}
				// Check whether or not individual is better than current best
				if (fit < bestfit)
				{
					bestfit = fit;
					bestid = id;
					for (int i = 0; i < NUM_TASKS; i++)
						bestTaskFit[i] = taskFit[id][i];
					for (int i = 0; i < glen; i++)
						bestInd[i] = m_pop[id][i];
				}
			}
            		exp->endGeneration(gn);
            		if ( commitStep() ) 
			{ 
				return; 
			}
			//printf("QUI\n");

			if (cstep < m_nevalsteps)
				reproduce(); // Reproduce
			else
			{
				// At the last generation we do not reproduce in order to keep the individuals of the last generation
				// Save the fitness and the number of visited cells by individuals of the last generation population
				sprintf(filename, "lastGenStatS%d.txt", currentSeed);
				fp = fopen(filename, "w");
				if (fp != NULL)
				{
					for(id = 0; id < popSize; id++)
						fprintf(fp, "%lf\n", tfitness[id] / ntfitness[id]);
					fclose(fp);
				}
				// Compute and save statistics
				computeFStat2();
				saveFStat();
				// Save best individual
				saveBestInd(bestInd);
				// Update generation counter
				cgen++;
			}

            		emit endGeneration( cgen, fmax, faverage, fmin );
            		if (commitStep()) 
			{
                		return; // stop the evolution process
            		}

            		farsa::Logger::info(QString("Generation %1 took %2 minutes - Generation's best fitness %3 - Best fitness = %4 - Evaluation steps = %5").arg(gn+1).arg((double)evotimer.elapsed()/60000.0, 0, 'f', 2).arg(fmin).arg(bestfit).arg(cstep));
            		fflush(stdout);
 
			//printf("QUI2\n");
			saveEvalStats(cstep, bestfit);

			// Save detailed stats about task fitnesses
			saveFitStats(bestTaskFit);
 
			//printf("QUI3\n");

			// Check whether the task is finished or not
			if (cstep >= m_nevalsteps)
			{
				if (!solved)
				{
					// Save the number of evaluation steps
					char evalfname[1024];
					sprintf(evalfname, "evalS%d.txt", currentSeed);
					FILE* evalfp = fopen(evalfname, "w");
					if (evalfp != NULL)
					{
						fprintf(evalfp, "%d", cstep);
						fclose(evalfp);
					}
					else
					{
						farsa::Logger::error("ERROR IN OPENING FILE " + QString(evalfname));
					}
				}
				
				tmpStatBuf += m_statBuf;
				// Saving files
				{
					const QString statFilename = QString("statS%1.fit").arg(currentSeed);
					QFile statFile(statFilename);
					if (!statFile.open(QIODevice::WriteOnly | QIODevice::Append | QIODevice::Truncate)) {
						farsa::Logger::warning("Cannot save stats into " + statFilename);
					} else {
						QTextStream s(&statFile);
						s << tmpStatBuf;
					}
					m_statBuf.clear();
				}

				{
					tmpFitBuf += m_fitBuf;
					// And best generalization fitness
					const QString fitFilename = QString("taskFitS%1.fit").arg(currentSeed);
					QFile fitFile(fitFilename);
					if (!fitFile.open(QIODevice::WriteOnly | QIODevice::Append | QIODevice::Truncate)) {
						farsa::Logger::warning("Cannot save stats into " + fitFilename);
					} else {
						QTextStream s(&fitFile);
						s << tmpFitBuf;
					}
					m_fitBuf.clear();
				}

				tmpEvalStatBuf += m_evalstatBuf;
				// Saving files
				{
					const QString evalStatFilename = QString("fitS%1.txt").arg(currentSeed);
					QFile evalStatFile(evalStatFilename);
					if (!evalStatFile.open(QIODevice::WriteOnly | QIODevice::Append | QIODevice::Truncate)) {
						farsa::Logger::warning("Cannot save stats into " + evalStatFilename);
					} else {
						QTextStream s(&evalStatFile);
						s << tmpEvalStatBuf;
					}
					//if (m_resetBufNumGen > 1)
					{
						m_evalstatBuf.clear();
					}
				}
			}
			
			// Saving files
			if (m_saveBestGenInfo)
			{
				if ((gn == 0) || (((gn + 1) % m_bestGenInterval) == 0))
				{
					const QString bestGenFilename = QString("B0S%1G%2.gen").arg(currentSeed).arg(gn + 1);
					QFile bestGenFile(bestGenFilename);
					if (!bestGenFile.open(QIODevice::WriteOnly | QIODevice::Append | QIODevice::Truncate)) {
						farsa::Logger::warning("Cannot save best genomes into " + bestGenFilename);
					} else {
						QTextStream s(&bestGenFile);
						s << m_bestGenBuf;
					}
					m_bestGenBuf.clear();
				}
			}
			gn++;
        	}
        	saveAllGen();
    	}
}

/* Sort samples based on fitness. Returns sorted fitnesses and corresponding
   indices. Flag <maximize> specifies the type of sorting (i.e., maximization or minimization) */
void MyEvoAlgo::sortSamples(QVector<double> f, const int size, QVector<double>& sf, QVector<int>& si, bool maximize)
{
	int i;
	double bf;
	int bi;
	int coeff;
	QVector<double> tmpf(size);
	if (maximize)
		coeff = -1;
	else
		coeff = 1;
	for (i = 0; i < size; i++)
		tmpf[i] = f[i];
	i = 0;
	while (i < size)
	{
		bf = tmpf[0];
		bi = 0;
		for (int j = 1; j < size; j++)
		{
			if (maximize)
			{
				// Maximization
				if (tmpf[j] > bf)
				{
					bf = tmpf[j];
					bi = j;
				}
			}
			else
			{
				// Minimization
				if (tmpf[j] < bf)
				{
					bf = tmpf[j];
					bi = j;
				}
			}
		}
		sf[i] = bf;
		si[i] = bi;
		tmpf[bi] = coeff * 999999999.0;
		i++;
	}
}

/*
 * Main function of the Genetic Algorithm (Steady State Version with the possibility
 * of performing the Annealing routine...)
 */
void MyEvoAlgo::evolveOpenAI()
{
	int rp; //replication
	int gn; //generation
	int cstep; // evaluation step
    	int id; //individuals
	double fit;
	int startGeneration = 0;
	char statfile[128];
	char genFile[128];
	char filename[64];
	FILE* fp;
    	int steps;
	bool solved;
	double bestfit;
	int bestid;
	QVector<float> bestTaskFit(NUM_TASKS);
	// Centroid
	QVector<float> center(glen);
	Eigen::VectorXf centroid(glen);
	Eigen::VectorXf dCen(glen);
	// Samples and offspring
	int lambda = 20;
	Eigen::MatrixXf samples(lambda, glen);
	Eigen::MatrixXf offspring(lambda * 2, glen);
	// Other variables
	Eigen::MatrixXf gradmat(1,glen);
	Eigen::VectorXf grad(glen); // Gradient
	Eigen::VectorXf globalg(glen); // Actual gradient!!!
	Eigen::VectorXf utilities(lambda * 2); // Array containing the utilities
	Eigen::VectorXf weights(lambda); // Array containing the weights
	Eigen::MatrixXf wmat(1,lambda);
	QVector<double> fitness(lambda * 2);
	QVector<double> sfitness(lambda * 2);
	QVector<int> sid(lambda * 2);
	// Vector to be used for evaluations
	QVector<float> ind(glen);
	// Noise std dev
	double noiseStdDev = 0.02;
	// Adam parameters
	double stepSize = 0.01; // Step size / learning rate
	double a;
	Eigen::VectorXf m(glen);
	Eigen::VectorXf v(glen);
	double beta1 = 0.9;
	double beta2 = 0.999;
	double epsilon = pow(10.0,-8.0);
	QVector<float> bestInd(glen);

	const int numInputs = (int)pow(2.0,(double)NUM_INPUTS);
	QVector<QVector<int> > bitStrings(numInputs);
	int** tmpStrings = generateBitstrings(NUM_INPUTS);
	for (int i = 0; i < numInputs; i++)
	{
		bitStrings[i].resize(NUM_INPUTS);
		for (int j = 0; j < NUM_INPUTS; j++)
		{
			bitStrings[i][j] = tmpStrings[i][j];
		}
	}

	QVector<float> taskFit(NUM_TASKS);
	
	farsa::Logger::info("EVOLUTION: OpenAI-ES with custom evaluation function");
    	farsa::Logger::info("Number of replications: " + QString::number(nreplications));

	// Individual to be evaluated
	QVector<int> evalInd(glen); // Vector to be used to avoid overwritings!!
	// replications
    	for(rp = 0; rp < nreplications; rp++) 
	{
		startGeneration = 0;
        	setSeed(getStartingSeed() + rp);
		farsa::Logger::info(QString("Replication %1, seed: %2").arg(rp + 1).arg(getStartingSeed() + rp));
		resetGenerationCounter();
		// Initialize the population (i.e. only the parents, offspring are mutated copies)
        	for (int j = 0; j < glen; j++)
		{
			center[j] = (float)getRng()->getDouble(-m_range, m_range);
			centroid[j] = center[j];
			m[j] = 0.0;
			v[j] = 0.0;
		}
		// Set fbest to a very low value
		this->fbest = 99999.0;

		emit startingReplication( rp );
       
		QTime evotimer;
        	
		// Buffers to save statistics, genotypes, fitnesses and the like
		m_bestGenBuf = "";
		m_statBuf = "";
		m_fitBuf = "";
		m_evalstatBuf = "";
		QString tmpStatBuf = "";
		QString tmpFitBuf = "";
		QString tmpEvalStatBuf = "";

		// Set the cart pole
		m_cartPole.setCartPole(getStartingSeed() + rp);
		// generations
		gn = startGeneration;
		cstep = 0;
		bestfit = 999999999.0;
		bestid = -1;
		solved = false;
		while (cstep < m_nevalsteps)
        	//for(gn=startGeneration;gn<nogenerations;gn++) 
		{
			evotimer.restart();
		    	farsa::Logger::info(" Generation " + QString::number(gn + 1));
		    	exp->initGeneration(gn);
		    	if ( commitStep() ) 
			{ 
				return;
			}

			// Create samples
			for (int i = 0; i < lambda; i++)
			{
				for (int j = 0; j < glen; j++)
				{
					samples(i,j) = getRng()->getGaussian(1.0, 0.0);
				}
			}
			// and offspring
			for (int i = 0; i < lambda; i++)
			{
				for (int j = 0; j < glen; j++)
				{
					for (int k = 0; k < 2; k++)
					{
						if (k % 2 == 0)
							offspring(i * 2 + k,j) = centroid[j] + noiseStdDev * samples(i,j);
						else
							offspring(i * 2 + k,j) = centroid[j] - noiseStdDev * samples(i,j);
					}
				}
			}
			
			//individuals
			for(int i = 0; i < lambda * 2; i++) 
			{
				// evaluate the parent
				//printf("### Offspring %d: ", i);
				for (int j = 0; j < glen; j++)
				{
					ind[j] = offspring(i,j);
					//printf("%f ", ind[j]);
				}
				steps = 0;
				fitness[i] = evaluate(ind, bitStrings, steps, taskFit);
				//printf("### - fitness = %lf\n", i, fitness[i]);
				if (isStopped()) 
				{ // stop evolution
					return;
				}
				cstep += steps;
				if ((fitness[i] == (double)NUM_TASKS) && !solved)
				{
					solved = true;
					// Save the number of evaluation steps
					char evalfname[1024];
					sprintf(evalfname, "evalS%d.txt", currentSeed);
					FILE* evalfp = fopen(evalfname, "w");
					if (evalfp != NULL)
					{
						fprintf(evalfp, "%d", cstep);
						fclose(evalfp);
					}
				}
				// Check whether or not individual is better than current best
				//printf("Sample %d: %lf vs %lf\n", i, fitness[i], bestfit);
				if (fitness[i] < bestfit)
				{
					bestfit = fitness[i];
					bestid = i;
					for (int j = 0; j < NUM_TASKS; j++)
						bestTaskFit[j] = taskFit[j];
					for (int j = 0; j < glen; j++)
						bestInd[j] = ind[j];
				}
			}
            		exp->endGeneration(gn);
            		if ( commitStep() ) 
			{ 
				return; 
			}

			cgen++;

			sortSamples(fitness, lambda * 2, sfitness, sid, false); // We need ascendent sorting!!!
			// Compute utilities
			for (int i = 0; i < lambda * 2; i++)
			{
				utilities[sid[i]] = (lambda * 2) - 1 - i; // Minimization: we need to invert utility scores
			}
			utilities /= ((lambda * 2) - 1);
			for (int i = 0; i < lambda * 2; i++)
			{
				utilities[i] -= 0.5;
			}
			// Compute weights
			for (int i = 0; i < lambda; i++)
			{
				weights[i] = (utilities[2 * i] - utilities[2 * i + 1]);
				wmat(0,i) = weights[i];
			}
			
			// Gradient
			gradmat = wmat * samples;

			for (int j = 0; j < glen; j++)
			{
				grad[j] = gradmat(0,j);
			}

			// Compute global gradient
			for (int j = 0; j < glen; j++)
			{
				globalg[j] = -grad[j];
				if (m_decay)
					globalg[j] += 0.005 * centroid[j];
			}
			
			// Adam optimization
			a = stepSize * sqrt(1.0 - pow(beta2, cgen)) / (1.0 - pow(beta1, cgen));
			// Update momentum vectors (mean and variance)
			for (int j = 0; j < glen; j++)
			{
				m[j] = beta1 * m[j] + (1.0 - beta1) * globalg[j];
				v[j] = beta2 * v[j] + (1.0 - beta2) * globalg[j] * globalg[j];
			}
			// Compute update of centroid
			for (int j = 0; j < glen; j++)
			{
				dCen[j] = -a * m[j] / (sqrt(v[j]) + epsilon);
				centroid[j] += dCen[j];
				center[j] = centroid[j];
			}

			if (m_refine)
			{
				QVector<float> original(glen);
				immediateMutation(center, original, 0);
				steps = 0;
				double cfit = evaluate(center, bitStrings, steps, taskFit);
				cstep += steps;
				QVector<float> novel(glen);
				QVector<float> novelTaskFit(NUM_TASKS);
				double novelFit = refine(original, cfit, taskFit, bitStrings, novel, cstep, novelTaskFit, solved);
				// Update the individual
				immediateMutation(novel, center, 0);
				for (int j = 0; j < glen; j++)
					centroid[j] = center[j];
				// Check whether or not the individual returned by the annealing routine is better than current best
				if (novelFit < bestfit)
				{
					bestfit = novelFit;
					bestid = -1;
					for (int i = 0; i < NUM_TASKS; i++)
						bestTaskFit[i] = novelTaskFit[i];
				}
			}

            		emit endGeneration( cgen, fmax, faverage, fmin );
            		if (commitStep()) 
			{
                		return; // stop the evolution process
            		}

            		farsa::Logger::info(QString("Generation %1 took %2 minutes - Generation's best fitness %3 - Best fitness = %4 - Evaluation steps = %5").arg(gn+1).arg((double)evotimer.elapsed()/60000.0, 0, 'f', 2).arg(sfitness[0]).arg(bestfit).arg(cstep));
            		fflush(stdout);
 
			//printf("QUI2\n");
			saveEvalStats(cstep, bestfit);

			// Save detailed stats about task fitnesses
			saveFitStats(bestTaskFit);
 
			//printf("QUI3\n");


			// Check whether the task is finished or not
			if (cstep >= m_nevalsteps)
			{
				if (!solved)
				{
					// Save the number of evaluation steps
					char evalfname[1024];
					sprintf(evalfname, "evalS%d.txt", currentSeed);
					FILE* evalfp = fopen(evalfname, "w");
					if (evalfp != NULL)
					{
						fprintf(evalfp, "%d", cstep);
						fclose(evalfp);
					}
					else
					{
						farsa::Logger::error("ERROR IN OPENING FILE " + QString(evalfname));
					}
				}
				
				tmpStatBuf += m_statBuf;
				// Saving files
				{
					const QString statFilename = QString("statS%1.fit").arg(currentSeed);
					QFile statFile(statFilename);
					if (!statFile.open(QIODevice::WriteOnly | QIODevice::Append | QIODevice::Truncate)) {
						farsa::Logger::warning("Cannot save stats into " + statFilename);
					} else {
						QTextStream s(&statFile);
						s << tmpStatBuf;
					}
					m_statBuf.clear();
				}

				{
					tmpFitBuf += m_fitBuf;
					// And best generalization fitness
					const QString fitFilename = QString("taskFitS%1.fit").arg(currentSeed);
					QFile fitFile(fitFilename);
					if (!fitFile.open(QIODevice::WriteOnly | QIODevice::Append | QIODevice::Truncate)) {
						farsa::Logger::warning("Cannot save stats into " + fitFilename);
					} else {
						QTextStream s(&fitFile);
						s << tmpFitBuf;
					}
					m_fitBuf.clear();
				}

				tmpEvalStatBuf += m_evalstatBuf;
				// Saving files
				{
					const QString evalStatFilename = QString("fitS%1.txt").arg(currentSeed);
					QFile evalStatFile(evalStatFilename);
					if (!evalStatFile.open(QIODevice::WriteOnly | QIODevice::Append | QIODevice::Truncate)) {
						farsa::Logger::warning("Cannot save stats into " + evalStatFilename);
					} else {
						QTextStream s(&evalStatFile);
						s << tmpEvalStatBuf;
					}
					//if (m_resetBufNumGen > 1)
					{
						m_evalstatBuf.clear();
					}
				}
			}
			
			// Saving files
			if (m_saveBestGenInfo)
			{
				if ((gn == 0) || (((gn + 1) % m_bestGenInterval) == 0))
				{
					const QString bestGenFilename = QString("B0S%1G%2.gen").arg(currentSeed).arg(gn + 1);
					QFile bestGenFile(bestGenFilename);
					if (!bestGenFile.open(QIODevice::WriteOnly | QIODevice::Append | QIODevice::Truncate)) {
						farsa::Logger::warning("Cannot save best genomes into " + bestGenFilename);
					} else {
						QTextStream s(&bestGenFile);
						s << m_bestGenBuf;
					}
					m_bestGenBuf.clear();
				}
			}
			gn++;
        	}
        	//saveAllGen();
		saveBestInd(bestInd);
    	}
}

/*
 * Main function of the Genetic Algorithm (Steady State Version with the possibility
 * of performing the Annealing routine...)
 */
void MyEvoAlgo::evolveNES()
{
	int rp; //replication
	int gn; //generation
	int cstep; // evaluation step
    	int id; //individuals
	double fit;
	int startGeneration = 0;
	char statfile[128];
	char genFile[128];
	char filename[64];
	FILE* fp;
    	int steps;
	bool solved;
	double bestfit;
	int bestid;
	QVector<float> bestTaskFit(NUM_TASKS);
	// Centroid
	QVector<float> center(glen);
	Eigen::VectorXf centroid(glen);
	Eigen::VectorXf dCen(glen);
	// Samples and offspring
	const float genLen = (const float)glen; // Floating point genome length
	float offspringSize = 4.0 + floor(3.0 * log(genLen)); // Floating point number of offspring generated by the individual
	const int offSize = (const int)floor(offspringSize); // Number of offspring generated by the individual
	Eigen::MatrixXf samples(offSize, glen);
	Eigen::MatrixXf offspring(offSize, glen);
	Eigen::MatrixXf trSamples(offSize, glen);
	Eigen::MatrixXf covMatrix(glen,glen); // Covariance Matrix
	Eigen::MatrixXf dCovM(glen,glen); // The update of the covariance matrix
	Eigen::MatrixXf expCovM(glen,glen); // Exponential matrix of CovM
	float etaInd = 1.0; // Learning rate for the individual
	float etaCovM = (9.0 + 3.0 * log(genLen)) / (5.0 * genLen * sqrt(genLen)); // Learning rate for the covariance matrix
	Eigen::VectorXf utility(offSize); // Array containing the utility values
	float utilSum; // Sum of utility values
	float tmp; // Temporary variable
	Eigen::VectorXf orderedUtility(offSize); // Array containing the ordered utilities (i.e. the best individual has the highest utility value)
	float ordUtilSum; // Sum of ordered utility values (it shall be equal to 1.0)
	Eigen::MatrixXf eye(glen,glen); // Identity matrix
	Eigen::MatrixXf ordUtil(glen,offSize);
	Eigen::MatrixXf ordUtilSamples(glen,offSize);
	Eigen::MatrixXf Z(offSize,glen);
	Eigen::MatrixXf G(glen,glen);
	Eigen::MatrixXf ordUtilTr(offSize,1);
	Eigen::MatrixXf v(glen,1);
	// Other variables
	QVector<double> fitness(offSize);
	QVector<double> sfitness(offSize);
	QVector<int> sid(offSize);
	// Vector to be used for evaluations
	QVector<float> ind(glen);
	QVector<float> bestInd(glen);

	//printf("glen %d genLen %f offspringSize %f offSize %d mu %d\n", glen, genLen, offspringSize, offSize, mu);
	
	utilSum = 0.0;
	// Calculate the utility values before normalization
	for (int j = 0; j < offSize; j++)
	{
		// u[j] is equal to the uttermost between 0 and (log(lambda / 2 + 1) - log(j))
		tmp = log(offspringSize / 2.0 + 1.0) - log((float)j + 1.0); // Remember that log(0) = -Inf, so 1 is added to index i
		if (tmp < 0.0)
		{
			tmp = 0.0;
		}
		utility[j] = tmp;
		// Update the sum of utility values
		utilSum += utility[j];
	}
	// Normalize the utility values
	for (int j = 0; j < offSize; j++)
	{
		if (utilSum != 0.0)
		{
			utility[j] /= utilSum;
		}
		utility[j] -= 1.0 / genLen; // In the standard version, utility values might be negative
	}

	for (int i = 0; i < glen; i++)
	{
		// Define the identity matrix
		eye(i,i) = 1.0;
		for (int j = 0; j < glen; j++)
		{
			if (i != j)
			{
				eye(i,j) = 0.0;
			}
			// Initialize covariance matrix, its update and its exponential
			covMatrix(i,j) = 0.0;
			dCovM(i,j) = 0.0;
			expCovM(i,j) = 0.0;
		}
	}

	const int numInputs = (int)pow(2.0,(double)NUM_INPUTS);
	QVector<QVector<int> > bitStrings(numInputs);
	int** tmpStrings = generateBitstrings(NUM_INPUTS);
	for (int i = 0; i < numInputs; i++)
	{
		bitStrings[i].resize(NUM_INPUTS);
		for (int j = 0; j < NUM_INPUTS; j++)
		{
			bitStrings[i][j] = tmpStrings[i][j];
		}
	}

	QVector<float> taskFit(NUM_TASKS);
	
	farsa::Logger::info("EVOLUTION: xNES with custom evaluation function");
    	farsa::Logger::info("Number of replications: " + QString::number(nreplications));

	// Individual to be evaluated
	QVector<int> evalInd(glen); // Vector to be used to avoid overwritings!!
	// replications
    	for(rp = 0; rp < nreplications; rp++) 
	{
		startGeneration = 0;
        	setSeed(getStartingSeed() + rp);
		farsa::Logger::info(QString("Replication %1, seed: %2").arg(rp + 1).arg(getStartingSeed() + rp));
		resetGenerationCounter();
		// Initialize the population (i.e. only the parents, offspring are mutated copies)
        	for (int j = 0; j < glen; j++)
		{
			center[j] = (float)getRng()->getDouble(-m_range, m_range);
			centroid[j] = center[j];
		}
		// Set fbest to a very low value
		this->fbest = 999999999.0;

		emit startingReplication( rp );
       
		QTime evotimer;
        	
		// Buffers to save statistics, genotypes, fitnesses and the like
		m_bestGenBuf = "";
		m_statBuf = "";
		m_fitBuf = "";
		m_evalstatBuf = "";
		QString tmpStatBuf = "";
		QString tmpFitBuf = "";
		QString tmpEvalStatBuf = "";

		// Initialize Covariance Matrix
		for (int i = 0; i < glen; i++)
		{
			for (int j = 0; j < glen; j++)
			{
				covMatrix(i,j) = 0.0;
				expCovM(i,j) = 0.0;
			}
		}

		// Set the cart pole
		m_cartPole.setCartPole(getStartingSeed() + rp);
		// generations
		gn = startGeneration;
		cstep = 0;
		bestfit = 999999999.0;
		bestid = -1;
		solved = false;
		while (cstep < m_nevalsteps)
        	//for(gn=startGeneration;gn<nogenerations;gn++) 
		{
			evotimer.restart();
		    	farsa::Logger::info(" Generation " + QString::number(gn + 1));
		    	exp->initGeneration(gn);
		    	if ( commitStep() ) 
			{ 
				return;
			}

			// Calculate the matrix exponential
			Eigen::MatrixExponential<Eigen::MatrixXf>(covMatrix).compute(expCovM);

			// Create samples
			for (int i = 0; i < offSize; i++)
			{
				for (int j = 0; j < glen; j++)
				{
					samples(i,j) = getRng()->getGaussian(1.0, 0.0);
				}
			}
			// and offspring
			trSamples = expCovM * samples;
			for (int i = 0; i < offSize; i++)
			{
				for (int j = 0; j < glen; j++)
				{
					offspring(i,j) = centroid[j] + trSamples(i,j);
				}
			}
			
			//individuals
			for(int i = 0; i < offSize; i++) 
			{
				// evaluate the parent
				//printf("### Offspring %d: ", i);
				for (int j = 0; j < glen; j++)
				{
					ind[j] = offspring(i,j);
					//printf("%f ", ind[j]);
				}
				steps = 0;
				fitness[i] = evaluate(ind, bitStrings, steps, taskFit);
				//printf("### %d - fitness = %lf\n", i, fitness[i]);
				if (isStopped()) 
				{ // stop evolution
					return;
				}
				cstep += steps;
				if ((fitness[i] == (double)NUM_TASKS) && !solved)
				{
					solved = true;
					// Save the number of evaluation steps
					char evalfname[1024];
					sprintf(evalfname, "evalS%d.txt", currentSeed);
					FILE* evalfp = fopen(evalfname, "w");
					if (evalfp != NULL)
					{
						fprintf(evalfp, "%d", cstep);
						fclose(evalfp);
					}
				}
				// Check whether or not individual is better than current best
				//printf("Sample %d: %lf vs %lf\n", i, fitness[i], bestfit);
				if (fitness[i] < bestfit)
				{
					bestfit = fitness[i];
					bestid = i;
					for (int j = 0; j < NUM_TASKS; j++)
						bestTaskFit[j] = taskFit[j];
					for (int j = 0; j < glen; j++)
						bestInd[j] = ind[j];
				}
			}
            		exp->endGeneration(gn);
            		if ( commitStep() ) 
			{ 
				return; 
			}

			cgen++;

			sortSamples(fitness, offSize, sfitness, sid, false); // We need ascendent sorting!!!

			// Initialize the sum of ordered utility values
			ordUtilSum = 0.0;
			// Now fill the array with the utility values
			for (int j = 0; j < offSize; j++)
			{
				orderedUtility[sid[j]] = utility[j];
				ordUtilSum += orderedUtility[sid[j]];
			}

			for (int i = 0; i < glen; i++)
			{
				for (int j = 0; j < offSize; j++)
				{
					ordUtil(i,j) = orderedUtility[j];
				}
			}

			for (int i = 0; i < glen; i++)
			{
				for (int j = 0; j < offSize; j++)
				{
					ordUtilSamples(i,j) = ordUtil(i,j) * samples(i,j);
					Z(j,i) = samples(i,j);
				}
			}

			G = ordUtilSamples * Z - ordUtilSum * eye;

			for (int j = 0; j < offSize; j++)
			{
				ordUtilTr(j,0) = orderedUtility[j];
			}

			v = etaInd * expCovM * samples * ordUtilTr;
			for (int i = 0; i < glen; i++)
			{
				dCen[i] = v(i,0);
			}
			dCovM = etaCovM * G;
			// Update the individual and the covariance matrix
			for (int j = 0; j < glen; j++)
			{
				centroid[j] += etaInd * dCen[j];
				center[j] = centroid[j];
			}
			covMatrix += dCovM;
			
            		emit endGeneration( cgen, fmax, faverage, fmin );
            		if (commitStep()) 
			{
                		return; // stop the evolution process
            		}

			/*double avgw = 0.0;
			for (int j = 0; j < glen; j++)
			{
				printf("%d: %lf\n", j, centroid[j]);
				avgw += centroid[j];
			}
			avgw /= glen;*/

            		farsa::Logger::info(QString("Generation %1 took %2 minutes - Generation's best fitness %3 - Best fitness = %4 - Evaluation steps = %5").arg(gn+1).arg((double)evotimer.elapsed()/60000.0, 0, 'f', 2).arg(sfitness[0]).arg(bestfit).arg(cstep));
            		fflush(stdout);
 
			//printf("QUI2\n");
			saveEvalStats(cstep, bestfit);

			// Save detailed stats about task fitnesses
			saveFitStats(bestTaskFit);
 
			//printf("QUI3\n");


			// Check whether the task is finished or not
			if (cstep >= m_nevalsteps)
			{
				if (!solved)
				{
					// Save the number of evaluation steps
					char evalfname[1024];
					sprintf(evalfname, "evalS%d.txt", currentSeed);
					FILE* evalfp = fopen(evalfname, "w");
					if (evalfp != NULL)
					{
						fprintf(evalfp, "%d", cstep);
						fclose(evalfp);
					}
					else
					{
						farsa::Logger::error("ERROR IN OPENING FILE " + QString(evalfname));
					}
				}
				
				tmpStatBuf += m_statBuf;
				// Saving files
				{
					const QString statFilename = QString("statS%1.fit").arg(currentSeed);
					QFile statFile(statFilename);
					if (!statFile.open(QIODevice::WriteOnly | QIODevice::Append | QIODevice::Truncate)) {
						farsa::Logger::warning("Cannot save stats into " + statFilename);
					} else {
						QTextStream s(&statFile);
						s << tmpStatBuf;
					}
					m_statBuf.clear();
				}

				{
					tmpFitBuf += m_fitBuf;
					// And best generalization fitness
					const QString fitFilename = QString("taskFitS%1.fit").arg(currentSeed);
					QFile fitFile(fitFilename);
					if (!fitFile.open(QIODevice::WriteOnly | QIODevice::Append | QIODevice::Truncate)) {
						farsa::Logger::warning("Cannot save stats into " + fitFilename);
					} else {
						QTextStream s(&fitFile);
						s << tmpFitBuf;
					}
					m_fitBuf.clear();
				}

				tmpEvalStatBuf += m_evalstatBuf;
				// Saving files
				{
					const QString evalStatFilename = QString("fitS%1.txt").arg(currentSeed);
					QFile evalStatFile(evalStatFilename);
					if (!evalStatFile.open(QIODevice::WriteOnly | QIODevice::Append | QIODevice::Truncate)) {
						farsa::Logger::warning("Cannot save stats into " + evalStatFilename);
					} else {
						QTextStream s(&evalStatFile);
						s << tmpEvalStatBuf;
					}
					//if (m_resetBufNumGen > 1)
					{
						m_evalstatBuf.clear();
					}
				}
			}
			
			// Saving files
			if (m_saveBestGenInfo)
			{
				if ((gn == 0) || (((gn + 1) % m_bestGenInterval) == 0))
				{
					const QString bestGenFilename = QString("B0S%1G%2.gen").arg(currentSeed).arg(gn + 1);
					QFile bestGenFile(bestGenFilename);
					if (!bestGenFile.open(QIODevice::WriteOnly | QIODevice::Append | QIODevice::Truncate)) {
						farsa::Logger::warning("Cannot save best genomes into " + bestGenFilename);
					} else {
						QTextStream s(&bestGenFile);
						s << m_bestGenBuf;
					}
					m_bestGenBuf.clear();
				}
			}
			gn++;
        	}
		saveBestInd(bestInd);
        	//saveAllGen();
    	}
}

/*
 * Main function of the Genetic Algorithm (Steady State Version with the possibility
 * of performing the Annealing routine...)
 */
void MyEvoAlgo::evolveNES2()
{
	int rp; //replication
	int gn; //generation
	int cstep; // evaluation step
    	int id; //individuals
	double fit;
	int startGeneration = 0;
	char statfile[128];
	char genFile[128];
	char filename[64];
	FILE* fp;
    	int steps;
	bool solved;
	double bestfit;
	int bestid;
	QVector<float> bestTaskFit(NUM_TASKS);
	// Centroid
	QVector<float> center(glen);
	Eigen::VectorXf centroid(glen);
	Eigen::VectorXf dCen(glen);
	// Samples and offspring
	const float genLen = (const float)glen; // Floating point genome length
	float offspringSize = 4.0 + floor(3.0 * log(genLen)); // Floating point number of offspring generated by the individual
	const int offSize = (const int)floor(offspringSize); // Number of offspring generated by the individual
	int mu = offSize / 2;
	Eigen::MatrixXf samples(offSize, glen);
	Eigen::MatrixXf offspring(offSize, glen);
	Eigen::MatrixXf sSamples(offSize, glen);
	Eigen::MatrixXf ssq(offSize,glen);
	Eigen::MatrixXf ones(offSize,glen);
	Eigen::MatrixXf ssqMinusOne(offSize,glen);
	// Other variables
	double initVar = 1.0;
	double centerLRate = 1.0;
	double covLRate = (3.0 + log(genLen)) / (5.0 * sqrt(genLen));
	double stepSize = 1.0 / (double)mu;
	Eigen::VectorXf weights(offSize);
	Eigen::MatrixXf wMat(1,offSize);
	Eigen::MatrixXf covGrad(1,glen);
	Eigen::MatrixXf dSigma(1,glen);
	Eigen::MatrixXf sigma(1,glen);
	QVector<double> fitness(offSize);
	QVector<double> sfitness(offSize);
	QVector<int> sid(offSize);
	// Vector to be used for evaluations
	QVector<float> ind(glen);
	QVector<float> bestInd(glen);

	//printf("glen %d genLen %f offspringSize %f offSize %d mu %d\n", glen, genLen, offspringSize, offSize, mu);
	
	// Init weights
	double wSum = 0.0;
	for (int i = 0; i < offSize; i++)
	{
		weights[i] = (log(offspringSize / 2.0 + 1.0) - log(i + 1));
		if (weights[i] < 0.0)
			weights[i] = 0.0;
		wSum += weights[i];
	}

	//double w = stepSize;
	// Set weights
	for (int i = 0; i < offSize; i++)
	{
		/*weights[offSize - mu + i] = w;
		wSum += weights[offSize - mu + i];
		w += stepSize;*/
		weights[i] /= wSum;
		weights[i] -= 1.0 / offspringSize;
		wMat(0,i) = weights[i];
		//printf("w%d %lf\n", i, weights[i]);
	}
	//exit(-1);

	/*for (int i = 0; i < offSize; i++)
	{
		weights[i] / wSum;
		wMat(0,i) = weights[i];
		//printf("w%d %lf\n", i, weights[i]);
	}*/

	// Init ones matrix
	for (int j = 0; j < glen; j++)
	{
		for (int i = 0; i < offSize; i++)
		{
			ones(i,j) = 1.0;
		}
	}

	const int numInputs = (int)pow(2.0,(double)NUM_INPUTS);
	QVector<QVector<int> > bitStrings(numInputs);
	int** tmpStrings = generateBitstrings(NUM_INPUTS);
	for (int i = 0; i < numInputs; i++)
	{
		bitStrings[i].resize(NUM_INPUTS);
		for (int j = 0; j < NUM_INPUTS; j++)
		{
			bitStrings[i][j] = tmpStrings[i][j];
		}
	}

	QVector<float> taskFit(NUM_TASKS);
	
	farsa::Logger::info("EVOLUTION: sNES with custom evaluation function");
    	farsa::Logger::info("Number of replications: " + QString::number(nreplications));

	// Individual to be evaluated
	QVector<int> evalInd(glen); // Vector to be used to avoid overwritings!!
	// replications
    	for(rp = 0; rp < nreplications; rp++) 
	{
		startGeneration = 0;
        	setSeed(getStartingSeed() + rp);
		farsa::Logger::info(QString("Replication %1, seed: %2").arg(rp + 1).arg(getStartingSeed() + rp));
		resetGenerationCounter();
		// Initialize the population (i.e. only the parents, offspring are mutated copies)
        	for (int j = 0; j < glen; j++)
		{
			center[j] = (float)getRng()->getDouble(-m_range, m_range);
			centroid[j] = center[j];
		}
		// Set fbest to a very low value
		this->fbest = 999999999.0;

		emit startingReplication( rp );
       
		QTime evotimer;
        	
		// Buffers to save statistics, genotypes, fitnesses and the like
		m_bestGenBuf = "";
		m_statBuf = "";
		m_fitBuf = "";
		m_evalstatBuf = "";
		QString tmpStatBuf = "";
		QString tmpFitBuf = "";
		QString tmpEvalStatBuf = "";

		// Init sigma
		for (int i = 0; i < glen; i++)
			sigma(0,i) = initVar;


		// Set the cart pole
		m_cartPole.setCartPole(getStartingSeed() + rp);
		// generations
		gn = startGeneration;
		cstep = 0;
		bestfit = 999999999.0;
		bestid = -1;
		solved = false;
		while (cstep < m_nevalsteps)
        	//for(gn=startGeneration;gn<nogenerations;gn++) 
		{
			evotimer.restart();
		    	farsa::Logger::info(" Generation " + QString::number(gn + 1));
		    	exp->initGeneration(gn);
		    	if ( commitStep() ) 
			{ 
				return;
			}

			// Create samples
			for (int i = 0; i < offSize; i++)
			{
				for (int j = 0; j < glen; j++)
				{
					samples(i,j) = getRng()->getGaussian(1.0, 0.0);
				}
			}
			// and offspring
			for (int i = 0; i < offSize; i++)
			{
				for (int j = 0; j < glen; j++)
				{
					offspring(i,j) = centroid[j] + 0.02 * sigma(0,j) * samples(i,j);
				}
			}
			
			//individuals
			for(int i = 0; i < offSize; i++) 
			{
				// evaluate the parent
				//printf("### Offspring %d: ", i);
				for (int j = 0; j < glen; j++)
				{
					ind[j] = offspring(i,j);
					//printf("%f ", ind[j]);
				}
				steps = 0;
				fitness[i] = evaluate(ind, bitStrings, steps, taskFit);
				//printf("### %d - fitness = %lf\n", i, fitness[i]);
				if (isStopped()) 
				{ // stop evolution
					return;
				}
				cstep += steps;
				if ((fitness[i] == (double)NUM_TASKS) && !solved)
				{
					solved = true;
					// Save the number of evaluation steps
					char evalfname[1024];
					sprintf(evalfname, "evalS%d.txt", currentSeed);
					FILE* evalfp = fopen(evalfname, "w");
					if (evalfp != NULL)
					{
						fprintf(evalfp, "%d", cstep);
						fclose(evalfp);
					}
				}
				// Check whether or not individual is better than current best
				//printf("Sample %d: %lf vs %lf\n", i, fitness[i], bestfit);
				if (fitness[i] < bestfit)
				{
					bestfit = fitness[i];
					bestid = i;
					for (int j = 0; j < NUM_TASKS; j++)
						bestTaskFit[j] = taskFit[j];
					for (int j = 0; j < glen; j++)
						bestInd[j] = ind[j];
				}
			}
            		exp->endGeneration(gn);
            		if ( commitStep() ) 
			{ 
				return; 

			}

			cgen++;

			sortSamples(fitness, offSize, sfitness, sid, false); // We need ascendent sorting!!!

			/*for (int i = 0; i < offSize; i++)
				printf("%d %d %lf\n", i, sid[i], sfitness[i]);
			exit(-1);*/

			for (int i = 0; i < offSize; i++)
			{
				int idx = sid[0];//offSize - i - 1];
				for (int j = 0; j < glen; j++)
				{
					sSamples(i,j) = samples(idx,j);
				}
			}

			dCen = wMat * sSamples;
			
			// Compute update of centroid
			for (int j = 0; j < glen; j++)
			{
				centroid[j] += centerLRate * dCen[j];
				center[j] = centroid[j];
			}

			// Update variances
			for (int i = 0; i < offSize; i++)
			{
				for (int j = 0; j < glen; j++)
				{
					ssq(i,j) = sSamples(i,j) * sSamples(i,j);
				}
			}
			ssqMinusOne = ssq - ones;
			covGrad = wMat * ssqMinusOne;
			dSigma = 0.5 * covLRate * covGrad;
			for (int j = 0; j < glen; j++)
			{
				sigma(0,j) *= std::exp(dSigma(0,j));
			}

            		emit endGeneration( cgen, fmax, faverage, fmin );
            		if (commitStep()) 
			{
                		return; // stop the evolution process
            		}


            		farsa::Logger::info(QString("Generation %1 took %2 minutes - Generation's best fitness %3 - Best fitness = %4 - Evaluation steps = %5").arg(gn+1).arg((double)evotimer.elapsed()/60000.0, 0, 'f', 2).arg(sfitness[0]).arg(bestfit).arg(cstep));
            		fflush(stdout);

			saveEvalStats(cstep, bestfit);
			saveFitStats(bestTaskFit);
			if (cstep >= m_nevalsteps)
			{
				if (!solved)
				{
					// Save the number of evaluation steps
					char evalfname[1024];
					sprintf(evalfname, "evalS%d.txt", currentSeed);
					FILE* evalfp = fopen(evalfname, "w");
					if (evalfp != NULL)
					{
						fprintf(evalfp, "%d", cstep);
						fclose(evalfp);
					}
					else
					{
						farsa::Logger::error("ERROR IN OPENING FILE " + QString(evalfname));
					}
				}

				tmpStatBuf += m_statBuf;
				// Saving files
				{
					const QString statFilename = QString("statS%1.fit").arg(currentSeed);
					QFile statFile(statFilename);
					if (!statFile.open(QIODevice::WriteOnly | QIODevice::Append | QIODevice::Truncate)) {
						farsa::Logger::warning("Cannot save stats into " + statFilename);
					} else {
						QTextStream s(&statFile);
						s << tmpStatBuf;
					}
					m_statBuf.clear();
				}

				{
					tmpFitBuf += m_fitBuf;
					// And best generalization fitness
					const QString fitFilename = QString("taskFitS%1.fit").arg(currentSeed);
					QFile fitFile(fitFilename);
					if (!fitFile.open(QIODevice::WriteOnly | QIODevice::Append | QIODevice::Truncate)) {
						farsa::Logger::warning("Cannot save stats into " + fitFilename);
					} else {
						QTextStream s(&fitFile);
						s << tmpFitBuf;
					}
					m_fitBuf.clear();
				}

				tmpEvalStatBuf += m_evalstatBuf;
				// Saving files
				{
					const QString evalStatFilename = QString("fitS%1.txt").arg(currentSeed);
					QFile evalStatFile(evalStatFilename);
					if (!evalStatFile.open(QIODevice::WriteOnly | QIODevice::Append | QIODevice::Truncate)) {
						farsa::Logger::warning("Cannot save stats into " + evalStatFilename);
					} else {
						QTextStream s(&evalStatFile);
						s << tmpEvalStatBuf;
					}
					//if (m_resetBufNumGen > 1)
					{
						m_evalstatBuf.clear();
					}
				}
			}
			
			// Saving files
			if (m_saveBestGenInfo)
			{
				if ((gn == 0) || (((gn + 1) % m_bestGenInterval) == 0))
				{
					const QString bestGenFilename = QString("B0S%1G%2.gen").arg(currentSeed).arg(gn + 1);
					QFile bestGenFile(bestGenFilename);
					if (!bestGenFile.open(QIODevice::WriteOnly | QIODevice::Append | QIODevice::Truncate)) {
						farsa::Logger::warning("Cannot save best genomes into " + bestGenFilename);
					} else {
						QTextStream s(&bestGenFile);
						s << m_bestGenBuf;
					}
					m_bestGenBuf.clear();
				}
			}
			gn++;
        	}
		saveBestInd(bestInd);
        	//saveAllGen();
    	}
}

void MyEvoAlgo::test()
{
	int rp; //replication
	int gn; //generation
	int id; //individuals
	QVector<double> bestFit(5);
	QVector<float> bestInd(glen);
	int startGeneration = 0;
	char statfile[128];
	char genFile[128];
	char filename[64];
	FILE* fp;

    	farsa::Logger::info("TEST: test evolved individuals on test functions");
    	farsa::Logger::info("Number of replications: " + QString::number(nreplications));

	// replications
	sprintf(filename, "testFunct.txt");
	fp = fopen(filename, "w");
	if (fp == NULL)
	{
		farsa::Logger::error(QString("Unable to open file ") + filename);
		exit(-1);
	}
    	for(rp = 0; rp < nreplications; rp++) 
	{
		startGeneration = 0;
        	setSeed(getStartingSeed() + rp);
		farsa::Logger::info(QString("Replication %1, seed: %2").arg(rp + 1).arg(getStartingSeed() + rp));
		resetGenerationCounter();

		// Initialize the population (i.e. only the parents, offspring are mutated copies)
        	initPop();

		loadBestInd(bestInd);

		getTestFitness(bestInd, bestFit);
		fprintf(fp, "%d\t", getStartingSeed() + rp);
		for (int i = 0; i < 5; i++)
			fprintf(fp, "%lf\t", bestFit[i]);
		fprintf(fp, "\n");
	}
	fclose(fp);
}

void MyEvoAlgo::getTestFitness(const QVector<float> ind, QVector<double>& fit)
{
	double fit_ra, fit_ro, fit_sp, fit_gr, fit_ac;
	if (fit.size() != 5)
		fit.resize(5);
	// Rastrigin
	fit_ra = rastrigin(ind);
	fit[0] = fit_ra;
	// Rosenbrock
	fit_ro = rosenbrock(ind);
	fit[1] = fit_ro;
	// Sphere
	fit_sp = sphere(ind);
	fit[2] = fit_sp;
	// Griewank
	fit_gr = griewank(ind);
    	fit[3] = fit_gr;
	// Ackley
	fit_ac = ackley(ind);
    	fit[4] = fit_ac;
}

int** MyEvoAlgo::generateBitstrings(int n){

    vector<string> input;
    int num_strings = (int)pow(2.0,(double)n);
    input.push_back("0");
    input.push_back("1");
    int i,j = 0;
    for(i = 2; i < (1<<n);i = i<<1){
        for(j = i-1; j>=0 ; j--)
            input.push_back(input[j]);

        for(j=0; j < i;j++)
            input[j] = "0" + input[j];

        for(j=i; j<2*i; j++)
            input[j] = "1" + input[j];
    }

    int** bit_input = new int*[num_strings];
    for(int i = 0; i< (int)pow(2.0,(double)n); i++)
        bit_input[i] = new int[n];


    for(i = 0; i< num_strings;i++){
        for(int j = 0; j < n; j++){
            string tmp ;
            tmp = input[i][j];
            bit_input[i][j] = atoi(tmp.c_str());
        }
    }
    return bit_input;

}

MyEvoAlgo::CartPole::CartPole()
	: MUP(0.000002)
	, MUC(0.0005)
	, GRAVITY(-9.8)
	, MASSCART(1.0)
	, MASSPOLE_1(0.1)
	, LENGTH_1(0.5)
	, FORCE_MAG(10.0)
	, TAU(0.01)
	, one_degree(0.0174532)
	, six_degrees(0.1047192)
	, twelve_degrees(0.2094384)
	, fifteen_degrees(0.2617993)
	, thirty_six_degrees(0.628329)
	, fifty_degrees(0.87266)
	, fixed_state(true)
	, rng(1)
{
}

MyEvoAlgo::CartPole::~CartPole()
{
}

void MyEvoAlgo::CartPole::initializeCartPole(bool velocity, bool fixedState){
  maxFitness = 100000;

  MARKOV=velocity;
  fixed_state = fixedState;
  stampOuput = "";

  LENGTH_2 = 0.05;
  MASSPOLE_2 = 0.01;

  // MyEvoAlgo::CartPole::reset() which is called here
}

void MyEvoAlgo::CartPole::setCartPole(const int s)
{
	// Set the seed
	setSeed(s);
}

farsa::RandomGenerator* MyEvoAlgo::CartPole::getRandomGenerator()
{
	return &rng;
}

void MyEvoAlgo::CartPole::setSeed(const int s)
{
	rng.setSeed(s);
}

double MyEvoAlgo::CartPole::evalNet(farsa::Evonet *evonet, int episode)
{
	int steps=0;
	double input[NUM_IN];
	double output;

	int nmarkovmax;

	double nmarkov_fitness;

	double jiggletotal; //total jiggle in last 100
	int count;  //step counter
	int ninputs; // number of inputs
	bool risk;

	if (nmarkov_long) nmarkovmax=100000;
	else nmarkovmax=1000;

	init(episode);

	if (MARKOV)
	{ // Velocities provided as inputs
		//printf("Episode %d - Initial state: %lf %lf %lf %lf %lf %lf\n", episode, state[0], state[1], state[2], state[3], state[4], state[5]);
		//while (steps++ < maxFitness)
		while (steps < nmarkovmax)
		{
			input[0] = state[0] / 4.8;
			input[1] = state[1] / 2;
			input[2] = state[2] / 0.52;
			input[3] = state[3] / 2;
			input[4] = state[4] / 0.52;
			input[5] = state[5] / 2;
			input[6] = 0.5;
			
			for (int j = 0; j < NUM_INPUTS; j++)
				evonet->setInput(j, input[j]);
			//Activate the net
			evonet->updateNet();
			// Get the output
			output = evonet->getOutput(0);
			// Update steps
			steps++;
			// Perform the action
			performAction(output,steps);
			// Check whether or not the pole(s) is(are) fallen
			if (outsideBounds())	// if failure
			{
				break;			// stop it now
			}
		}
		return ((double) steps);
	}
	else
	{  //NON MARKOV CASE
		while (steps < nmarkovmax)
		{
			input[0] = state[0] / 4.8;
			input[1] = state[2] / 0.52;
			input[2] = state[4] / 0.52;
			//input[3] = 0.5;
			risk = alert();
			if (risk)
				input[3] = 1.0;
			else
				input[3] = 0.0;
			ninputs = 4;
			
			for (int j = 0; j < ninputs; j++)
				evonet->setInput(j, input[j]);
			//Activate the net
			evonet->updateNet();
			// Get the output
			output=evonet->getOutput(0);
			// Update steps
			steps++;
			// Perform the action
			performAction(output,steps);
			// Check whether or not the pole(s) is(are) fallen
			if (outsideBounds())	// if failure
				break;			// stop it now

			if (nmarkov_long && (outsideBounds()))	// if failure
				break;			// stop it now
		}
		return (double) steps;
	}
}

void MyEvoAlgo::CartPole::init(int episode)
{
	int j;

	if (!MARKOV) {
		//Clear all fitness records
		cartpos_sum=0.0;
		cartv_sum=0.0;
		polepos_sum=0.0;
		polev_sum=0.0;
	}

	balanced_sum=0; //Always count # balanced

	last_hundred=false;

	// Evaluation phase
	if (fixed_state)
	{
		// Fixed initial states
		if (episode < 0)
		{
			state[0] = state[1] = state[3] = state[4] = state[5] = 0;
			state[2] = 0.07; // Approximately 4°
			//farsa::Logger::info(QString("state: %1 %2 %3 %4 0.5").arg(state[0]).arg(state[1]).arg(state[2]).arg(state[3]));
		}
		else
		{
			if (episode == 0)
			{
				state[0] = -1.944;
				state[1] = state[2] = state[3] = state[4] = state[5] = 0;
			}
			if (episode == 1)
			{
				state[0] = 1.944;
				state[1] = state[2] = state[3] = state[4] = state[5] = 0;
			}
			if (episode == 2)
			{
				state[1] = -1.215;
				state[0] = state[2] = state[3] = state[4] = state[5] = 0;
			}
			if (episode == 3)
			{
				state[1] = 1.215;
				state[0] = state[2] = state[3] = state[4] = state[5] = 0;
			}
			if (episode == 4)
			{
				state[2] = -0.10472;
				state[0] = state[1] = state[3] = state[4] = state[5] = 0;
			}
			if (episode == 5)
			{
				state[2] = 0.10472;
				state[0] = state[1] = state[3] = state[4] = state[5] = 0;
			}
			if (episode == 6)
			{
				state[3] = -0.135088;
				state[0] = state[1] = state[2] = state[4] = state[5] = 0;
			}
			if (episode == 7)
			{
				state[3] = 0.135088;
				state[0] = state[1] = state[2] = state[4] = state[5] = 0;
			}
		}
		//Sleep(100000);
	}
	else
	{
		// Random initial states
		state[0] = getRandomGenerator()->getDouble(-1.944,1.944);
		state[1] = getRandomGenerator()->getDouble(-1.215,1.215);
		state[2] = getRandomGenerator()->getDouble(-0.10472,0.10472);
		state[3] = getRandomGenerator()->getDouble(-0.135088,0.135088);
		state[4] = getRandomGenerator()->getDouble(-0.10472,0.10472);
		state[5] = getRandomGenerator()->getDouble(-0.135088,0.135088);
	}
}

double MyEvoAlgo::CartPole::getTaskLength()
{
	if (nmarkov_long)
		return 100000.0;
	else
		return 1000.0;
}

void MyEvoAlgo::CartPole::performAction(double output, int stepnum)
{

	int i;
	double  dydx[6];

	const bool RK4 = true; //Set to Runge-Kutta 4th order integration method
	const double EULER_TAU = TAU / 4;
  
	/*--- Apply action to the simulated cart-pole ---*/
	if(RK4)
	{
		for(i=0;i<2;++i)
		{
			dydx[0] = state[1];
			dydx[2] = state[3];
			dydx[4] = state[5];
			step(output,state,dydx);
			rk4(output,state,dydx,state);
		}
	}
	else
	{
		for(i=0;i<8;++i)
		{
			dydx[0] = state[1];
			dydx[2] = state[3];
			dydx[4] = state[5];
			step(output,state,dydx);
			state[0] += EULER_TAU * dydx[0];
			state[1] += EULER_TAU * dydx[1];
			state[2] += EULER_TAU * dydx[2];
			state[3] += EULER_TAU * dydx[3];
			state[4] += EULER_TAU * dydx[4];
			state[5] += EULER_TAU * dydx[5];
		}
	}

	//Record this state
	cartpos_sum += fabs(state[0]);
	cartv_sum += fabs(state[1]);
	polepos_sum += fabs(state[2]);
	polev_sum += fabs(state[3]);
	if (stepnum <= 1000)
		jigglestep[stepnum-1]=fabs(state[0])+fabs(state[1])+fabs(state[2])+fabs(state[3]);

	if (!(outsideBounds()))
		++balanced_sum;
}

void MyEvoAlgo::CartPole::step(double action, double *st, double *derivs)
{
    	double force,costheta_1,costheta_2,sintheta_1,sintheta_2,
          gsintheta_1,gsintheta_2,temp_1,temp_2,ml_1,ml_2,fi_1,fi_2,mi_1,mi_2;
	
	force =  (action - 0.5) * FORCE_MAG * 2;
	costheta_1 = cos(st[2]);
	sintheta_1 = sin(st[2]);
	gsintheta_1 = GRAVITY * sintheta_1;
	costheta_2 = cos(st[4]);
	sintheta_2 = sin(st[4]);
	gsintheta_2 = GRAVITY * sintheta_2;

	ml_1 = LENGTH_1 * MASSPOLE_1;
	ml_2 = LENGTH_2 * MASSPOLE_2;
	temp_1 = MUP * st[3] / ml_1;
	temp_2 = MUP * st[5] / ml_2;
	fi_1 = (ml_1 * st[3] * st[3] * sintheta_1) +
			(0.75 * MASSPOLE_1 * costheta_1 * (temp_1 + gsintheta_1));
	fi_2 = (ml_2 * st[5] * st[5] * sintheta_2) +
			(0.75 * MASSPOLE_2 * costheta_2 * (temp_2 + gsintheta_2));
	mi_1 = MASSPOLE_1 * (1 - (0.75 * costheta_1 * costheta_1));
	mi_2 = MASSPOLE_2 * (1 - (0.75 * costheta_2 * costheta_2));

	derivs[1] = (force + fi_1 + fi_2)
					/ (mi_1 + mi_2 + MASSCART);

	derivs[3] = -0.75 * (derivs[1] * costheta_1 + gsintheta_1 + temp_1)
					/ LENGTH_1;
	derivs[5] = -0.75 * (derivs[1] * costheta_2 + gsintheta_2 + temp_2)
					/ LENGTH_2;
}

void MyEvoAlgo::CartPole::rk4(double f, double y[], double dydx[], double yout[])
{

    int i;

    double hh,h6,dym[6],dyt[6],yt[6];


    hh=TAU*0.5;
    h6=TAU/6.0;
    for (i=0;i<=5;i++) yt[i]=y[i]+hh*dydx[i];
    step(f,yt,dyt);
    dyt[0] = yt[1];
    dyt[2] = yt[3];
    dyt[4] = yt[5];
    for (i=0;i<=5;i++) yt[i]=y[i]+hh*dyt[i];
    step(f,yt,dym);
    dym[0] = yt[1];
    dym[2] = yt[3];
    dym[4] = yt[5];
    for (i=0;i<=5;i++) {
        yt[i]=y[i]+TAU*dym[i];
        dym[i] += dyt[i];
    }
    step(f,yt,dyt);
    dyt[0] = yt[1];
    dyt[2] = yt[3];
    dyt[4] = yt[5];
    for (i=0;i<=5;i++)
        yout[i]=y[i]+h6*(dydx[i]+dyt[i]+2.0*dym[i]);
}

bool MyEvoAlgo::CartPole::outsideBounds()
{
	const double failureAngle = thirty_six_degrees;

	return
		(state[0] < -2.4             ||
		state[0] > 2.4               ||
		state[2] < -failureAngle     ||
		state[2] > failureAngle      ||
		state[4] < -failureAngle     ||
		state[4] > failureAngle);
}

bool MyEvoAlgo::CartPole::alert()
{
	return (fabs(state[0]) >= 2.0 || fabs(state[2]) >= 0.5 || fabs(state[4]) >= 0.5);
}

MyEvoAlgo::Grid::Grid()
	: size(GRID_SIZE)
{
}

MyEvoAlgo::Grid::~Grid()
{
}

double MyEvoAlgo::Grid::evalNet(farsa::Evonet *evonet, int trial, double& dist)
{
	double input[4];
	int steps = 0;
	double output;
	bool collision = false;
	bool reached = false;

	init(trial);

	//printf("Trial %d\n", trial);

	while (steps < GRID_STEPS)
	{
		// Convert states into inputs in [0,1]
		/*input[0] = state[0] / (GRID_SIZE - 1) - 0.5;
		input[1] = state[1] / (GRID_SIZE - 1) - 0.5;
		input[2] = 1.0 - calcDistFromTarget() / GRID_SIZE;//distFromTarget() / sqrt(GRID_SIZE * GRID_SIZE);
		input[3] = angleFromTarget() / PI_GRECO;*/
                for (int j = 0; j < 4; j++)
			input[j] = state[j] / (GRID_SIZE - 1); // Values in the range [-0.5,0.5]
		evonet->resetNet();
		for (int j = 0; j < 4; j++)
			evonet->setInput(j, input[j]);
		//Activate the net
		evonet->updateNet();
		// Get the output
		output = evonet->getOutput(0);
		// Perform the action and check for collisions with walls
		collision = performAction(output);
		// Check whether the target has been reached
		reached = reachTarget();
		// Update steps
		steps++;
		// Stop episode in case of collisions or if the target has been reached
		if (collision || reached)
			break;
	}
	
	dist = calcDistFromTarget();
	
	return steps;
}

void MyEvoAlgo::Grid::init(int trial)
{
	// Target position is the center of the grid
	state[2] = (double)((GRID_SIZE - 1) / 2);
	state[3] = (double)((GRID_SIZE - 1) / 2);

	// Agent position depends on the episode
	if (trial == 0)
	{
		// Bottom-left corner
		state[0] = 0.0;
		state[1] = 0.0;
	}
	else if (trial == 1)
	{
		// Top-left corner
		state[0] = 0.0;
		state[1] = (double)(GRID_SIZE - 1);
	}
	else if (trial == 2)
	{
		// Top-right corner
		state[0] = (double)(GRID_SIZE - 1);
		state[1] = (double)(GRID_SIZE - 1);
	}
	else
	{
		// Bottom-right corner
		state[0] = (double)(GRID_SIZE - 1);
		state[1] = 0.0;
	}
	//printf("Starting agent position: (%lf,%lf)\n", state[0], state[1]);
}

bool MyEvoAlgo::Grid::performAction(double output)
{
	bool collision = false;
	int dx, dy;
	int i;

	if ((output >= -1.0) && (output < -0.5))
	{
		// Left
		dx = -1;
		dy = 0;
	}
	else if ((output >= -0.5) && (output < 0.0))
	{
		// Top
		dx = 0;
		dy = 1;
	}
	else if ((output >= 0.0) && (output < 0.5))
	{
		// Right
		dx = 1;
		dy = 0;
	}
	else
	{
		// Bottom
		dx = 0;
		dy = -1;
	}
	
	// Update agent position and check boundaries
	state[0] += (double)dx;
	state[1] += (double)dy;
	if ((state[0] < 0.0) || (state[0] >= (double)GRID_SIZE))
	{
		if (state[0] < 0.0)
			state[0] = 0.0;
		if (state[0] >= (double)GRID_SIZE)
			state[0] = (double)(GRID_SIZE - 1);
		collision = true;
	}
	if ((state[1] < 0.0) || (state[1] >= (double)GRID_SIZE))
	{
		if (state[1] < 0.0)
			state[1] = 0.0;
		if (state[1] >= (double)GRID_SIZE)
			state[1] = (double)(GRID_SIZE - 1);
		collision = true;
	}
	//printf("New agent position: (%lf,%lf)\n", state[0], state[1]);
	return collision;
}

bool MyEvoAlgo::Grid::reachTarget()
{
	bool reached = false;
	//printf("agent pos: (%lf,%lf), target pos: (%lf,%lf)\n", state[0], state[1], state[2], state[3]);
	if ((state[0] == state[2]) && (state[1] == state[3]))
		reached = true;
	return reached;
}

double MyEvoAlgo::Grid::calcDistFromTarget()
{
	double dx = fabs(state[2] - state[0]);
	double dy = fabs(state[3] - state[1]);
	return dx + dy;
}

double MyEvoAlgo::Grid::distFromTarget()
{
	double dx = (state[0] - state[2]);
	double dy = (state[1] - state[3]);
	double dist = sqrt((dx * dx) + (dy * dy));
	return dist;
}

double MyEvoAlgo::Grid::angleFromTarget()
{
	double dx = (state[2] - state[0]);
	double dy = (state[3] - state[1]);
	if (dx == 0.0 && dy == 0.0)
		return 0.0;
	return atan2(dy, dx);
}

MyEvoAlgo::Arena::Arena()
	: size(ARENA_SIZE)
{
}

MyEvoAlgo::Arena::~Arena()
{
}

double MyEvoAlgo::Arena::evalNet(farsa::Evonet *evonet, int& num)
{
	double input[4];
	int steps = 0;
	double output;
	bool collision = false;

	init();

	//printf("Trial %d\n", trial);

	while (steps < ARENA_STEPS)
	{
		// Convert states into inputs
		input[0] = state[0] / (ARENA_SIZE - 1) - 0.5;
		input[1] = state[1] / (ARENA_SIZE - 1) - 0.5;
		input[2] = state[2];
		input[3] = state[3];
		evonet->resetNet();
		for (int j = 0; j < 4; j++)
			evonet->setInput(j, input[j]);
		//Activate the net
		evonet->updateNet();
		// Get the output
		output = evonet->getOutput(0);
		// Perform the action and check for collisions with walls
		collision = performAction(output);
		// Update steps
		steps++;
		// Stop episode in case of collisions or if the target has been reached
		if (collision)
			break;
	}
	
	// Calculate the number of visited cells
	num = calcMat();
	
	return steps;
}

void MyEvoAlgo::Arena::init()
{
	int x, y, num;

	// Initialize matrix
	initMat();
	
	// Agent position is the center of the arena
	x = (ARENA_SIZE - 1) / 2;
	y = (ARENA_SIZE - 1) / 2;
	state[0] = (double)x;
	state[1] = (double)y;

	// Third state is the ground sensor
	state[2] = mat[x][y];

	// Fourth state is the number of already visited cells
	num = 0;
	state[3] = (double)num / (ARENA_SIZE * ARENA_SIZE);

	// Update current cell to visited
	mat[x][y] = 1.0;
}

void MyEvoAlgo::Arena::initMat()
{
	for (int i = 0; i < ARENA_SIZE; i++)
	{
		for (int j = 0; j < ARENA_SIZE; j++)
		{
			mat[i][j] = 0.0;
		}
	}
}

int MyEvoAlgo::Arena::calcMat()
{
	int num = 0;
	for (int i = 0; i < ARENA_SIZE; i++)
	{
		for (int j = 0; j < ARENA_SIZE; j++)
		{
			if (mat[i][j] == 1.0)
				num++;
		}
	}
	return num;
}

bool MyEvoAlgo::Arena::performAction(double output)
{
	bool collision = false;
	int x, y, dx, dy;
	int num;

	if ((output >= -1.0) && (output < -0.5))
	{
		// Left
		dx = -1;
		dy = 0;
	}
	else if ((output >= -0.5) && (output < 0.0))
	{
		// Top
		dx = 0;
		dy = 1;
	}
	else if ((output >= 0.0) && (output < 0.5))
	{
		// Right
		dx = 1;
		dy = 0;
	}
	else
	{
		// Bottom
		dx = 0;
		dy = -1;
	}
	
	// Update agent position and check boundaries
	state[0] += (double)dx;
	state[1] += (double)dy;
	if ((state[0] < 0.0) || (state[0] >= (double)ARENA_SIZE))
	{
		if (state[0] < 0.0)
			state[0] = 0.0;
		if (state[0] >= (double)ARENA_SIZE)
			state[0] = (double)(ARENA_SIZE - 1);
		collision = true;
	}
	if ((state[1] < 0.0) || (state[1] >= (double)ARENA_SIZE))
	{
		if (state[1] < 0.0)
			state[1] = 0.0;
		if (state[1] >= (double)ARENA_SIZE)
			state[1] = (double)(ARENA_SIZE - 1);
		collision = true;
	}
	// Get agent pos
	x = (int)state[0];
	y = (int)state[1];
	// Update ground sensor
	state[2] = mat[x][y];
	// Update number of visited cells
	num = calcMat();
	state[3] = (double)num / (ARENA_SIZE * ARENA_SIZE);
	// Update matrix
	if (mat[x][y] == 0.0)
		mat[x][y] = 1.0;
	//printf("New agent position: (%lf,%lf)\n", state[0], state[1]);
	return collision;
}

MyEvoAlgo::Game::Game()
	: width(GAME_WIDTH)
	, height(GAME_HEIGHT)
	, rng(1)
{
}

MyEvoAlgo::Game::~Game()
{
}

double MyEvoAlgo::Game::evalNet(farsa::Evonet *evonet, int trial)
{
	double input[4];
	int steps = 0;
	double output;
	double score = 0.0;

	init(trial);

	while (steps < GAME_STEPS)
	{
		// Convert states into inputs in [0,1]
		input[0] = (state[0] / (GAME_WIDTH - 1) - 0.5) * 2.0;
		input[1] = (state[1] / (GAME_HEIGHT - 1) - 0.5) * 2.0;
		input[2] = (state[2] / (GAME_WIDTH - 1) - 0.5) * 2.0;
		input[3] = (state[3] / (GAME_HEIGHT - 1) - 0.5) * 2.0;
		evonet->resetNet();
		for (int j = 0; j < 4; j++)
			evonet->setInput(j, input[j]);
		//Activate the net
		evonet->updateNet();
		// Get the output
		output = evonet->getOutput(0);
		// Perform the action
		performAction(output);
		// Update steps
		steps++;
		// Check whether the agent avoid object falling
		if (state[3] == 0.0)
		{
			//printf("Check state: %lf %lf\n", state[0], state[2]);
			if (state[0] == state[2])
				// Agent reached the object before it falls
				score = 1.0;
			break;
		}
	}
	return score;
}

void MyEvoAlgo::Game::init(int trial)
{
	// x-coord of the falling object is random
	state[2] = (double)trial;
	// y-coord of the falling object is the height of the game
	state[3] = (double)(GAME_HEIGHT - 1);
	// x-coord of the agent is random
	state[0] = (double)((GAME_WIDTH - 1) / 2.0);
	// y-coord of the agent is 0
	state[1] = 0.0;
	//printf("State: %lf %lf %lf %lf\n", state[0], state[1], state[2], state[3]);
}

void MyEvoAlgo::Game::performAction(double output)
{
	// Object falls
	state[3] -= 1.0;

	// Agent can move either left or right
	if (output <= 0.0)
		// Left
		state[0] -= 1.0;
	else
		// Right
		state[0] += 1.0;
	
	// Check boundaries
	if (state[0] < 0.0)
		state[0] = 0.0;
	if (state[0] > (double)(GAME_WIDTH - 1))
		state[0] = (double)(GAME_WIDTH - 1);
}

farsa::RandomGenerator* MyEvoAlgo::Game::getRandomGenerator()
{
	return &rng;
}

void MyEvoAlgo::Game::setSeed(const int s)
{
	rng.setSeed(s);
}

