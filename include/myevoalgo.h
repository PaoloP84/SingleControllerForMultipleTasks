/********************************************************************************
 *  FARSA - Total99                                                             *
 *  Copyright (C) 2012-2013 Gianluca Massera <emmegian@yahoo.it>                *
 *                                                                              *
 *  This program is free software; you can redistribute it and/or modify        *
 *  it under the terms of the GNU General Public License as published by        *
 *  the Free Software Foundation; either version 2 of the License, or           *
 *  (at your option) any later version.                                         *
 *                                                                              *
 *  This program is distributed in the hope that it will be useful,             *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of              *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the               *
 *  GNU General Public License for more details.                                *
 *                                                                              *
 *  You should have received a copy of the GNU General Public License           *
 *  along with this program; if not, write to the Free Software                 *
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA  *
 ********************************************************************************/

#ifndef MYEVOALGO_H
#define MYEVOALGO_H

#include "farsaplugin.h"

#include "evoga.h"
#include "logger.h"
#include "evodataviewer.h"
#include "randomgenerator.h"
#include "evorobotexperiment.h"
#include "cartpoleexperiment.h"
#include "evonet.h"
#include <parametersettable.h>
#include <configurationhelper.h>
#include <configurationparameters.h>

#include <vector>
#include <fstream>
#include <iostream>
#include <deque>
#include <iomanip>
#include <sstream>
#include <string>

#include <QFile>
#include <QTime>
#include <QVector>
#include <QtAlgorithms>
#include <QThreadPool>
#include <QtConcurrentMap>

#define NUM_TASKS 4 // 4-bit parity, double-pole balancing, grid navigation, function approximation, game
#define NUM_INPUTS 4 // 4 inputs for all tasks
#define DPOLE_EPS 8 // Fixed Initial States condition (see Pagliuca et al. (2018))
#define GRID_SIZE 501
#define GRID_EPS 4 // From each corner
#define GRID_STEPS 1000
#define FUNCT_EPS 12 // Cases with pseudo-random numbers in the range [-5.0,5.0]
#define ARENA_SIZE 21
#define ARENA_STEPS 500
#define GAME_WIDTH 21
#define GAME_HEIGHT 50
#define GAME_EPS GAME_WIDTH
#define GAME_STEPS GAME_HEIGHT
#define LOCAL_SEED 9966

/**
 * \brief My evolutionary algorithm. It inherits from Evoga
 *        and runs a modified version of the steadyState
 *        algorithm (with floating-point genes).
 */
class FARSA_PLUGIN_API MyEvoAlgo : public farsa::Evoga
{
	Q_OBJECT
	FARSA_REGISTER_CLASS(Evoga)

public:
     /**
	 * \brief Constructor
	 */
	 MyEvoAlgo();


     /**
	 * \brief Destructor
	 */     
    ~MyEvoAlgo();
	
	/**
	 * \brief Configures the object using a ConfigurationParameters object
	 *
	 * \param params the configuration parameters object with parameters to
	 *               use
	 * \param prefix the prefix to use to access the object configuration
	 *               parameters. This is guaranteed to end with the
	 *               separator character when called by the factory, so you
	 *               don't need to add one
	 */
	virtual void configure(farsa::ConfigurationParameters& params, QString prefix);

	virtual void save(farsa::ConfigurationParameters& params, QString prefix)
	{
		// Calling parent function
		farsa::Evoga::save(params, prefix);
		farsa::Logger::warning("NOT IMPLEMENTED (MyEvoAlgo::save)");
		abort();
	}

	//! Describe
	static void describe(QString type);
	
	void postConfigureInitialization();

	/**
	 * \brief It invokes either evolve() or test() method independently of
	 *        the chosen type. This depends on the <test> flag.
	 */
	virtual void evolveAllReplicas();

	void setSeed(const int seed);

private:
	/**
	 * \brief Entry point of MyEvoAlgo class
	 */
	void evolveSSS();
	/**
	 * \brief Entry point of MyEvoAlgo class
	 */
	void evolveGGA();
	/**
	 * \brief Entry point of MyEvoAlgo class
	 */
	void evolveOpenAI();
	/**
	 * \brief Entry point of MyEvoAlgo class
	 */
	void evolveNES();
	/**
	 * \brief Entry point of MyEvoAlgo class
	 */
	void evolveNES2();
	/**
	 * \brief Entry point of test method (functions only!!!)
	 */
	void test();
	/**
	 * \brief Initializes the population
	 */
	void initPop();
	/**
	 * \brief Resets the population (new size = 0)
	 */
	void resetPop();
	/**
	 * \brief Mutates the gene
	 *
	 * \param from the source individual
	 * \param to the destination individual
	 * \param mut the mutation rate
	 */
	void performMutation(int from, int to, int mut);
	/**
	 * \brief Mutates the gene (does not consider the population)
	 *
	 * \param from the source individual
	 * \param to the destination individual
	 * \param mut the mutation rate
	 */
	void immediateMutation(const QVector<float> from, QVector<float>& to, int mut);
	/**
	 * \brief Neutral muttaion (does not consider the population)
	 *
	 * \param from the source individual
	 * \param to the destination individual
	 */
	void neutralMutation(const QVector<float> from, QVector<float>& to);
	/**
	 * \brief Perform the evaluation of the individual
	 *        and return the corresponding fitness
	 *
	 * \param ind the individual to be evaluated
	 * \return the fitness of the individual
	 */
	double evaluateParity(const QVector<float> ind, const QVector< QVector<int> > bitStrings, int& steps);
	/**
	 * \brief Perform the evaluation of the individual
	 *        and return the corresponding fitness
	 *
	 * \param ind the individual to be evaluated
	 * \return the fitness of the individual
	 */
	double evaluateDpole(const QVector<float> ind, int& steps);
	/**
	 * \brief Perform the evaluation of the individual
	 *        and return the corresponding fitness
	 *
	 * \param ind the individual to be evaluated
	 * \return the fitness of the individual
	 */
	double evaluateGrid(const QVector<float> ind, int& steps);
	/**
	 * \brief Perform the evaluation of the individual
	 *        and return the corresponding fitness
	 *
	 * \param ind the individual to be evaluated
	 * \return the fitness of the individual
	 */
	double evaluateFunct(const QVector<float> ind, int& steps);
	/**
	 * \brief Perform the evaluation of the individual
	 *        and return the corresponding fitness
	 *
	 * \param ind the individual to be evaluated
	 * \return the fitness of the individual
	 */
	double evaluateExploration(const QVector<float> ind, int& steps);
	/**
	 * \brief Perform the evaluation of the individual
	 *        and return the corresponding fitness
	 *
	 * \param ind the individual to be evaluated
	 * \return the fitness of the individual
	 */
	double evaluateGame(const QVector<float> ind, int& steps);
	/**
	 * \brief Perform the evaluation of the individual
	 *        and return the corresponding fitness
	 *
	 * \param ind the individual to be evaluated
	 * \return the fitness of the individual
	 */
	double evaluate(const QVector<float> ind, const QVector< QVector<int> > bitStrings, int& steps, QVector<float>& taskFit);
	/**
	 * \brief Loads a genotype
	 *
	 * \param fp the file descriptor
	 * \param ind the individual to be loaded
	 */
	void loadGen(FILE *fp, int ind);
	/**
	 * \brief Loads all the genotypes
	 *
	 * \param gen the generation
	 * \param filew the filename
	 * \return the number of loaded individuals
	 */
	int loadAllGen(int gen, char* filew);
	/**
	 * \brief Saves the best individual

	 *
	 */
	void loadBestInd(QVector<float>& ind);
	/**
	 * \brief Saves the best individual
	 *
	 */
	void saveBestInd(QVector<float> ind);
	/**
	 * \brief Saves a genotype
	 *
	 * \param fp The file pointer
	 * \param ind The individual to be saved
	 */
	void saveGen(FILE *fp, int ind);
	/**
	 * \brief Saves all the genotypes
	 *
	 */
	void saveAllGen();
	/**
	 * \brief Computes fitness statistics (min, max, avg)
	 *
	 */
	void computeFStat2();
	/**
	 * \brief Saves fitness statistics (min, max, avg)
	 *
	 */
	void saveFStat();
	/**
	 * \brief Saves detailed fitness stats (for each task) and the steps
	 *
	 */
	void saveFitStats(QVector<float> fit);
	/**
	 * \brief Saves evaluation stats
	 *
	 */
	void saveEvalStats(int steps, double fit);
	/**
	 * \brief Fills a QVector with -1 values (invalid)
	 *
	 * \param ind The QVector
	 */
	void invalidFill(QVector<float>& ind);
	/**
	 * \brief Refinement routine (it could be useless!!!)
	 *
	 * \param original The original individual (before the refinement routine)
	 * \param novel The new individual (after the refinement routine)
	 * \return The fitness of the novel individual
	 */
	double refine(const QVector<float> original, double fitness, QVector<float> taskFit, const QVector< QVector<int> > bitStrings, QVector<float>& novel, int& steps, QVector<float>& novelTaskFit, bool& solved);

	//! Reproduction with mutations and asexual crossover
	void reproduce();

	//! Sort samples
	void sortSamples(QVector<double> f, const int size, QVector<double>& sf, QVector<int>& si, bool maximize);

	//! Max values of optimization test functions
	void calcMaxVals();
	//! Rastrigin function
	double rastrigin(QVector<float> x);
	//! Rosenbrock function
	double rosenbrock(QVector<float> x);
	//! Sphere function
	double sphere(QVector<float> x);
	//! Griewank function
	double griewank(QVector<float> x);
	//! Ackley function
	double ackley(QVector<float> x);

	//! Get fitness values
	void getTestFitness(const QVector<float> ind, QVector<double>& fit);

	//! Compute input patterns of length <n>
	int** generateBitstrings(int n);

	//! Get random number generator
	farsa::RandomGenerator* getRng();

	//! Get random number generator
	farsa::RandomGenerator* getRngFunct();

	//! Experiment
	farsa::CartPoleExperiment* m_exp;
	//! Population
	QVector< QVector<float> > m_pop;
	//! Number of evaluation steps
	int m_nevalsteps;
	//! Mutation rate
	double m_mutationRate;
	//! Crossover rate
	double m_crossoverRate;

	//! Weight range
	double m_range;

	//! Algorithm id
	int m_algo; // 0: Stochastic Steady State (SSS), 1: Generational Genetic Algorithm (GGA) with elitism
	bool m_decay;

	//! Best selected (for GGA)
	QVector< QVector<float> > m_bests;
	
	//! Buffers for saving best genotypes and best fitness
	QString m_bestGenBuf;
	QString m_statBuf;
	QString m_fitBuf;
	QString m_evalstatBuf;
	
	//! Flags whether or not the information about the best individual must be saved
	bool m_saveBestGenInfo;
	//! Interval for saving information about the best individual
	int m_bestGenInterval;
	//! Noise on the individual's fitness (i.e., stochasticity)
	double m_fitNoise;
	//! Flags whether or not the refinement must be run
	bool m_refine;
	//! Number of refinement iterations
	int m_numRefineIterations;
	//! Flags whether or not randomization must be performed in order to avoid the replacement of the same worst individuals
	//! Random generator
	farsa::RandomGenerator m_rng;
	//! Random generator for function approximation task
	farsa::RandomGenerator m_rngFunct;

	//! Test functions only
	bool m_test;
	
	//! The class implementing the Cart-Pole
	class CartPole {
	public:
		CartPole();
		~CartPole();
		void initializeCartPole(bool velocity, bool fixedState);
		void setCartPole(const int s);
		virtual double evalNet(farsa::Evonet* evonet, int thresh);
		/**
		 * \brief Returns the random number generator
		 *
		 * \return The random number generator
		 */
		farsa::RandomGenerator* getRandomGenerator();
		void setSeed(const int s);
		double getTaskLength();
		double maxFitness;
		bool MARKOV;
		bool fixed_state;
		
		//! Local random number generator
		farsa::RandomGenerator rng;
		
		bool last_hundred;
		bool nmarkov_long;  //Flag that we are looking at the champ
		QString stampOuput;
		
		double state[6];
		
		double jigglestep[1000];
	
	protected:
		virtual void init(int type);
	
	private:
	
		void performAction(double output, int stepnum);
		void step(double action, double *state, double *derivs);
		void rk4(double f, double y[], double dydx[], double yout[]);
		bool outsideBounds();
		bool alert();
		
		const static int NUM_IN=7;
		//const static double MUP = 0.000002;
		//const static double MUC = 0.0005;
		//const static double GRAVITY= -9.8;
		//const static double MASSCART= 1.0;
		//const static double MASSPOLE_1= 0.1;
		
		//const static double LENGTH_1= 0.5;		  /* actually half the pole's length */
		
		//const static double FORCE_MAG= 10.0;
		//const static double TAU= 0.01;		  //seconds between state updates
		
		//const static double one_degree= 0.0174532;	/* 2pi/360 */
		//const static double six_degrees= 0.1047192;
		//const static double twelve_degrees= 0.2094384;
		//const static double fifteen_degrees= 0.2617993;
		//const static double thirty_six_degrees= 0.628329;
		//const static double fifty_degrees= 0.87266;

		const double MUP;
		const double MUC;
		const double GRAVITY;
		const double MASSCART;
		const double MASSPOLE_1;
		
		const double LENGTH_1;		  /* actually half the pole's length */
		
		const double FORCE_MAG;
		const double TAU;		  //seconds between state updates
		
		const double one_degree;	/* 2pi/360 */
		const double six_degrees;
		const double twelve_degrees;
		const double fifteen_degrees;
		const double thirty_six_degrees;
		const double fifty_degrees;
		
		double LENGTH_2;
		double MASSPOLE_2;
		
		//Queues used for Gruau's fitness which damps oscillations
		int balanced_sum;
		double cartpos_sum;
		double cartv_sum;
		double polepos_sum;
		double polev_sum;
	};
	
	// The cart pole used for evaluating individuals
	CartPole m_cartPole;

	class Grid {
	public:
		Grid();
		~Grid();
		void init(int trial);
		virtual double evalNet(farsa::Evonet* evonet, int trial, double& dist);

		// Grid size
		int size;
		
		double state[4];
	
	private:
	
		bool performAction(double output);
		bool reachTarget();
		double distFromTarget();
		double angleFromTarget();
		double calcDistFromTarget();
	};

	// Grid environment
	Grid m_grid;

	class Arena {
	public:
		Arena();
		~Arena();
		void init();
		virtual double evalNet(farsa::Evonet* evonet, int& num);

		// Arena size
		int size;
		
		double state[4];

		// Matrix storing cells
		double mat[ARENA_SIZE][ARENA_SIZE];
	
	private:
	
		bool performAction(double output);
		void initMat();
		int calcMat();
	};

	// Arena environment
	Arena m_arena;

	class Game {
	public:
		Game();
		~Game();
		void init(int trial);
		virtual double evalNet(farsa::Evonet* evonet, int trial);
		/**
		 * \brief Returns the random number generator
		 *
		 * \return The random number generator
		 */
		farsa::RandomGenerator* getRandomGenerator();
		void setSeed(const int s);

		//! Local random number generator
		farsa::RandomGenerator rng;
		
		// Environment width and height
		int width;
		int height;
		
		double state[4];
	
	private:
	
		void performAction(double output);
	};

	//! Game environment
	Game m_game;
};

#endif
