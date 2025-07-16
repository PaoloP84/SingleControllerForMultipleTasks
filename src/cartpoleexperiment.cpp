/********************************************************************************
 *  FARSA - Total99                                                             *
 *  Copyright (C) 2005-2011 Gianluca Massera <emmegian@yahoo.it>                *
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

#include "cartpoleexperiment.h"
#include "utilitiesexceptions.h"
#include "randomgenerator.h"
#include "robots.h"
#include "world.h"
#include "phybox.h"
#include "wheeledexperimenthelper.h"
#include "logger.h"
#include "configurationhelper.h"
#include "helperresources.h"

// This is needed because the isnan and isinf functions are not present windows
#ifdef WIN32
	#include <float.h>
	#define isnan(x) _isnan(x)
	#define isinf(x) (!_finite(x))
#else
	#define isnan(x) std::isnan(x)
	#define isinf(x) std::isinf(x)
#endif

#define LOCAL_SEED 7531

namespace farsa
{

CartPoleExperiment::CartPoleExperiment() :
	farsa::EvoRobotExperiment(),
	m_recreateWorld(false),
	m_robots(),
	m_robotNames(),
	m_playgroundCellSize(0.1f),
    	m_playgroundCellMatrix(),
    	m_cellMap(),
    	m_cellMapSize(0),
	m_numVisitedCells(0),
	m_rng(LOCAL_SEED)
{
	// No need to add resources, they are added by EvoRobotExperiment
	addUsableResource("agent[0]:FakeInput");
}

CartPoleExperiment::~CartPoleExperiment()
{
}

void CartPoleExperiment::configure(farsa::ConfigurationParameters& params, QString prefix)
{
	// Calling parent function
	farsa::EvoRobotExperiment::configure(params, prefix);
	m_playgroundCellSize = farsa::ConfigurationHelper::getDouble(params, prefix + "playgroundCellSize", m_playgroundCellSize);
}

void CartPoleExperiment::save(farsa::ConfigurationParameters& params, QString prefix)
{
	// Calling parent function
	farsa::EvoRobotExperiment::save(params, prefix);

    farsa::Logger::error("NOT IMPLEMENTED (CartPoleExperiment::save)");
	abort();
}

void CartPoleExperiment::describe(QString type)
{
	// Calling parent function
	farsa::EvoRobotExperiment::describe(type);

	Descriptor d = addTypeDescription(type, "The cart-pole experiment");
}

void CartPoleExperiment::postConfigureInitialization()
{
	// Calling parent function
	EvoRobotExperiment::postConfigureInitialization();

	farsa::ResourcesLocker locker(this);
	farsa::Arena* arena = getResource<farsa::Arena>("arena");
	double size = ((arena->getWidth() * arena->getHeight()) / (m_playgroundCellSize * m_playgroundCellSize));

    	m_cellMapSize = ceil(size);

    	// Find out the size of the matrix (i.e. rows and columns)
    	const int rows = (int)floor(arena->getWidth() / m_playgroundCellSize);
    	const int cols = (int)floor(arena->getHeight() / m_playgroundCellSize);
    	m_playgroundCellMatrix.resize(rows,cols);

    	// Initialize cell map (i.e. the map that associates a cell to an interval)
    	initCellMap();

	// Here we create the arena and the object
    	setupArena();
}

void CartPoleExperiment::initGeneration(int /*generation*/)
{
}

void CartPoleExperiment::initIndividual(int /*individual*/)
{
	// Resetting the world to avoid numerical problems
	m_recreateWorld = true;

	setSeed(LOCAL_SEED); // Deterministic task
}

void CartPoleExperiment::initTrial(int /*trial*/)
{
	if (m_recreateWorld) {
        	recreateWorld();
        	recreateAllRobots();
        	recreateArena();

        	setupArena();

        	m_recreateWorld = false;
    	}

	farsa::ResourcesLocker locker(this);
    	farsa::Arena* arena = getResource<farsa::Arena>("arena");

	farsa::RobotOnPlane* robot = getResource<farsa::RobotOnPlane>("agent[0]:robot");

	// Robot is placed here!!!
	const farsa::real arenaWidthHalfLimitForRobot = arena->getWidth() / 2.0 - robot->robotRadius() - 0.1;
	const farsa::real arenaHeightHalfLimitForRobot = arena->getHeight() / 2.0 - robot->robotRadius() - 0.1;
	const farsa::real rx = getRng()->getDouble(-arenaWidthHalfLimitForRobot, arenaWidthHalfLimitForRobot);
	const farsa::real ry = getRng()->getDouble(-arenaHeightHalfLimitForRobot, arenaHeightHalfLimitForRobot);
	const farsa::real ro = getRng()->getDouble(-PI_GRECO, PI_GRECO);
	robot->setPosition(arena->getPlane(), rx, ry);
	robot->setOrientation(arena->getPlane(), ro);

	// Resetting fitness and event counters for the current trial
    	trialFitnessValue = 0.0;

	// Reset number of visited cells
	m_numVisitedCells = 0;

	// Initialise the matrix
	for (int i = 0; i < m_playgroundCellMatrix.rows(); i++)
    	{
		for (int j = 0; j < m_playgroundCellMatrix.cols(); j++)
		{
			m_playgroundCellMatrix(i,j) = 0.0;
		}
    	}

	// Reset neural network
	farsa::Evonet* evonet = getResource<farsa::Evonet>("agent[0]:evonet");
	evonet->resetNet();
}

void CartPoleExperiment::initStep(int /*step*/)
{
	farsa::ResourcesLocker locker(this);
	farsa::RobotOnPlane* robot = getResource<farsa::RobotOnPlane>("agent[0]:robot");
	// Set robot sensors
	setSensors(robot);
}

void CartPoleExperiment::afterSensorsUpdate()
{
}

void CartPoleExperiment::beforeMotorsUpdate()
{
}

void CartPoleExperiment::beforeWorldAdvance()
{
	farsa::ResourcesLocker locker(this);
	farsa::RobotOnPlane* robot = getResource<farsa::RobotOnPlane>("agent[0]:robot");
	// Get robot speed limits
	double minSp1, minSp2, maxSp1, maxSp2;
	farsa::PhyEpuck* r = dynamic_cast<farsa::PhyEpuck*>(robot);
	// Get network output
	farsa::Evonet* evonet = getResource<farsa::Evonet>("agent[0]:evonet");
	double output = evonet->getOutput(0);
	// Convert the output into motor speeds
	QVector<double> v(2);
	if (output < 0.25)
		// Move forward
		v[0] = v[1] = 1.0;
	else if (output >= 0.25 && output < 0.5)
	{
		// Turn left
		v[0] = 0.0;
		v[1] = 1.0;
	}
	else if (output >= 0.5 && output < 0.75)
	{
		// Move backward
		v[0] = -1.0;
		v[1] = -1.0;
	}
	else
	{
		// Turn right
		v[0] = 1.0;
		v[1] = 0.0;
	}
	r->wheelsController()->setSpeeds(maxSp1 * v[0], maxSp2 * v[1]);
}

void CartPoleExperiment::endStep(int /*step*/)
{
	farsa::ResourcesLocker locker(this);
    	farsa::RobotOnPlane* robot = getResource<farsa::RobotOnPlane>("agent[0]:robot");
	farsa::Arena* arena = getResource<farsa::Arena>("arena");

	// Check collisions with either food objects or walls
	QVector< farsa::PhyObject2DWrapper * > collidedObjects = arena->getKinematicRobotCollisions("agent[0]:robot");
	int n = collidedObjects.size();
	if (n > 0)
	{
		for (int i = 0; i < n; i++) 
		{
			if ( collidedObjects[i]->type() == farsa::PhyObject2DWrapper::Wall ) 
			{
				// Set random orientation
				robot->setOrientation(arena->getPlane(), getRng()->getDouble(-PI_GRECO, PI_GRECO));
			}
		}
	}

	// Update number of visited cells
	int cellIdx = getCellIndex(robot);
	const int x = m_cellMap[cellIdx].rowIdx;
	const int y = m_cellMap[cellIdx].colIdx;
	if (m_playgroundCellMatrix(y,x) != 1.0)
	{
		m_playgroundCellMatrix(y,x) = 1.0;
	}
}

void CartPoleExperiment::endTrial(int /*trial*/)
{
	// Count the number of visited cells
	m_numVisitedCells = 0;
	for (int i = 0; i < m_playgroundCellMatrix.rows(); i++)
	{
		for (int j = 0; j < m_playgroundCellMatrix.cols(); j++)
		{
			if (m_playgroundCellMatrix(i,j) == 1.0)
				m_numVisitedCells++;
		}
	}
	trialFitnessValue = (double)m_numVisitedCells / (m_playgroundCellMatrix.rows() * m_playgroundCellMatrix.cols());
	totalFitnessValue += trialFitnessValue;
}

void CartPoleExperiment::endIndividual(int /*individual*/)
{
	totalFitnessValue = totalFitnessValue / farsa::real(getNTrials());
}

void CartPoleExperiment::endGeneration(int /*generation*/)
{
}

void CartPoleExperiment::resourceChanged(QString resourceName, ResourceChangeType changeType)
{
	// Calling parent function
	farsa::EvoRobotExperiment::resourceChanged(resourceName, changeType);

	// Here we are only interested in robots, so we build a regular expression to only check when
	// a robot changes
	QRegExp checkRobots("agent\\[(\\d)\\]:robot");
	if (checkRobots.indexIn(resourceName) != -1) {
		// Getting the index
		int index = checkRobots.cap(1).toUInt();

		// The code below should work as expected. If there is a crash, then there is probably a bug
		// in resource management
		if (changeType == Deleted) {
			// Removing robot
			m_robots[index] = NULL;
		} else {
			// Adding the robot to our list
			farsa::RobotOnPlane* const robot = getResource<farsa::RobotOnPlane>();

			if (m_robots.size() == index) {
				m_robots.append(robot);
				m_robotNames.append(QString("agent["+QString::number(index)+"]:robot"));
			} else {
				m_robots[index] = robot;
				m_robotNames[index] = QString("agent["+QString::number(index)+"]:robot");
			}
		}
	}
}

void CartPoleExperiment::setupArena()
{
    	farsa::ResourcesLocker locker(this);

    	// Getting the arena
    	farsa::Arena* arena = getResource<farsa::Arena>("arena");

	//arena->getPlane()->setColor(Qt::black);

    	// Creating walls all around the arena
    	const real m_wallThickness = 0.01f;
    	const real m_wallHeight = 0.05f;
    	const real halfHeight = arena->getHeight() / 2.0;
    	const real halfWidth = arena->getWidth() / 2.0;
    	const real topWallPos = halfHeight + m_wallThickness / 2.0;
    	const real bottomWallPos = -(halfHeight + m_wallThickness / 2.0);
    	const real rightWallPos = halfWidth + m_wallThickness / 2.0;
    	const real leftWallPos = -(halfWidth + m_wallThickness / 2.0);
    	arena->createWall(Qt::black, wVector(-halfWidth, topWallPos, 0.0), wVector(halfWidth, topWallPos, 0.0), m_wallThickness, m_wallHeight);
    	arena->createWall(Qt::black, wVector(-halfWidth, bottomWallPos, 0.0), wVector(halfWidth, bottomWallPos, 0.0), m_wallThickness, m_wallHeight);
    	arena->createWall(Qt::black, wVector(rightWallPos, -halfHeight, 0.0), wVector(rightWallPos, halfHeight, 0.0), m_wallThickness, m_wallHeight);
    	arena->createWall(Qt::black, wVector(leftWallPos, -halfHeight, 0.0), wVector(leftWallPos, halfHeight, 0.0), m_wallThickness, m_wallHeight);

    	// Finally checking if the arena is big enough. We only check against the robot radius because it is bigger than the object radius
    	farsa::RobotOnPlane* robot = getResource<farsa::RobotOnPlane>("agent[0]:robot");
    	const real minDimension = robot->robotRadius() * 10;
    	if (arena->getWidth() < minDimension) 
	{
        	farsa::throwUserRuntimeError("Cannot run the experiment because the arena is too small (width too small)");
    	} 
	else if (arena->getHeight() < minDimension) 
	{
        	farsa::throwUserRuntimeError("Cannot run the experiment because the arena is too small (height too small)");
    	}
}

void CartPoleExperiment::setSensors(farsa::RobotOnPlane* robot)
{
	farsa::ResourcesLocker locker(this);
	farsa::Arena* arena = getResource<farsa::Arena>("arena");
	ResourceVector<farsa::real>* fakeInput = getResource< ResourceVector<farsa::real> > ("agent[0]:FakeInput");
	const int size = fakeInput->size();
	if (size != 4)
	{
		farsa::Logger::error(QString("Invalid size %1 for energy sensor... Abort!!!").arg(size));
		exit(-1);
	}
	const farsa::real halfHeight = arena->getHeight() / 2.0;
	const farsa::real halfWidth = arena->getWidth() / 2.0;
	//! Position and orientation of the robot
	farsa::wVector robotPos = robot->position();
	double robotAng = robot->orientation(arena->getPlane()) + PI_GRECO / 2.0; // To account for differences between orientation() and setOrientation() methods
	// Normalize angle in the range [-pi,pi]
	if (robotAng < -PI_GRECO)
		robotAng += 2.0 * PI_GRECO;
	if (robotAng > PI_GRECO)
		robotAng -= 2.0 * PI_GRECO;
	//! Set sensors
	//! Normalized position and orientation (sensors in the range [0,1])
	(*fakeInput)[0] = (1.0 + robotPos.x / halfWidth) / 2.0;
	(*fakeInput)[1] = (1.0 + robotPos.y / halfHeight) / 2.0;
	(*fakeInput)[2] = (1.0 + robotAng / PI_GRECO) / 2.0;
	//! Last sensor works as a ground sensor!!!
	int cellIdx = getCellIndex(robot);
	const int x = m_cellMap[cellIdx].rowIdx;
	const int y = m_cellMap[cellIdx].colIdx;
	if (m_playgroundCellMatrix(y,x) == 1.0)
		(*fakeInput)[3] = 1.0;
	else
		(*fakeInput)[3] = 0.0;
}

int CartPoleExperiment::getCellIndex(farsa::RobotOnPlane* robot)
{
	// Getting the robot position
	const farsa::wVector robotPos = robot->position();
	// To find out the cell indices we use the information stored in the cell map
	bool found = false;
	int cellIdx = -1; // Invalid value!!!
	int i = 0;
	// First, we look for the x index of the cell
	while ((i < m_cellMapSize) && !found)
	{
		// Get the ranges for x and y
		float minx = m_cellMap[i].bottomleftx;
		float maxx = m_cellMap[i].uprightx;
		float miny = m_cellMap[i].bottomlefty;
		float maxy = m_cellMap[i].uprighty;
		//printf("min x y %lf %lf - max x y %lf %lf\n", minx, miny, maxx, maxy);
		// Check if the robot is inside the current cell
		if (((robotPos.x >= minx) && (robotPos.x < maxx)) && ((robotPos.y >= miny) && (robotPos.y < maxy)))
		{
			cellIdx = i;
			found = true;
		}
		i++;
	}
	// If the robot is not in any of cells, stop the execution!
	if (!found)
	{
		farsa::Logger::info("Robot position .x: " + QString::number(robotPos.x) + ", .y: " + QString::number(robotPos.y));
		farsa::Logger::error("FATAL ERROR!!! cell index not found! Check your code...");
		exit(-1);
	}

	return cellIdx;
}

void CartPoleExperiment::initCellMap()
{
	farsa::ResourcesLocker locker(this);
    	farsa::Arena* arena = getResource<farsa::Arena>("arena");

	// Getting the sizes of the arena
	const farsa::real halfHeight = arena->getHeight() / 2.0;
	const farsa::real halfWidth = arena->getWidth() / 2.0;
	float x = halfWidth;
	float y = halfHeight;
	const int rows = (int)floor(arena->getWidth() / m_playgroundCellSize);
	const int cols = (int)floor(arena->getHeight() / m_playgroundCellSize);
	int rowIdx = 0;
	int colIdx = 0;
	// Compute the size of the map
	m_cellMap.resize(m_cellMapSize);
	int i;
	// Initialise the map (to avoid numerical issues)
	for (i = 0; i < m_cellMapSize; i++)
	{
		m_cellMap[i].rowIdx = 0;
		m_cellMap[i].colIdx = 0;
		m_cellMap[i].bottomleftx = 0.0;
		m_cellMap[i].bottomlefty = 0.0;
		m_cellMap[i].uprightx = 0.0;
		m_cellMap[i].uprighty = 0.0;
	}
	// Fill the map
	i = 0;
	while (i < m_cellMapSize)
	{
		m_cellMap[i].rowIdx = rowIdx;
		m_cellMap[i].colIdx = colIdx;
		m_cellMap[i].bottomleftx = x - m_playgroundCellSize;
		m_cellMap[i].bottomlefty = y - m_playgroundCellSize;
		m_cellMap[i].uprightx = x;
		m_cellMap[i].uprighty = y;
		// Update
		x -= m_playgroundCellSize;
		colIdx++;
		i++;
		if ((i % rows) == 0)
		{
			// Switch row
			rowIdx++;
			colIdx = 0;
			x = halfWidth;
			y -= m_playgroundCellSize;
		}
	}
}

cell2interval CartPoleExperiment::mapCellToInterval(int rowIdx, int colIdx)
{
	// Given row and column indices, this function returns the pointer to the struct
	// containing information about bottom-left and up-right vertices
	cell2interval retVal;
	int i = 0;
	bool found = false;
	// Compute the size of the map
	// Find the entry of the map containing the required information
	while ((i < m_cellMapSize) && !found)
	{
		if ((rowIdx == m_cellMap[i].rowIdx) && (colIdx == m_cellMap[i].colIdx))
		{
			retVal = m_cellMap[i];
			found = true;
		}
		i++;
	}
	// If none of the elements of the map contains the required information, stop the execution!
	if (!found)
	{
		farsa::Logger::error(QString("Unable to find the element of the map corresponding to the given indices %1 %2!").arg(rowIdx).arg(colIdx));
		exit(-1);
	}
	return retVal;
}

void CartPoleExperiment::createCells()
{
	farsa::ResourcesLocker locker(this);
	farsa::Arena* arena = getResource<farsa::Arena>("arena");
	// Find out the size of the matrix (i.e. rows and columns)
	const int rows = (int)floor(arena->getWidth() / m_playgroundCellSize);
	const int cols = (int)floor(arena->getHeight() / m_playgroundCellSize);
	m_playgroundCells.resize(rows * cols);
	for (int j = 0; j < cols; j++)
	{
		for (int i = 0; i < rows; i++)
		{
			// Create a "rectangular" (it is a square) area
			m_playgroundCells[(i * cols) + j] = arena->createRectangularTargetArea(m_playgroundCellSize, m_playgroundCellSize, QColor(0,0,0,0));
			cell2interval mapEntry = mapCellToInterval(j, i);
			const float x = (mapEntry.bottomleftx + mapEntry.uprightx) / 2.0;
			const float y = (mapEntry.bottomlefty + mapEntry.uprighty) / 2.0;
			// Set the position of the rectangular area
			m_playgroundCells[(i * cols) + j]->setPosition(x,y);
		}
	}
}

void CartPoleExperiment::destroyCells()
{
	farsa::ResourcesLocker locker(this);
	farsa::Arena* arena = getResource<farsa::Arena>("arena");
	// Find out the size of the matrix (i.e. rows and columns)
	const int rows = (int)floor(arena->getWidth() / m_playgroundCellSize);
	const int cols = (int)floor(arena->getHeight() / m_playgroundCellSize);
	// Deleting all cells, they will be added again below
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			arena->delete2DObject(m_playgroundCells[(i * rows) + j]);
		}
	}
	m_playgroundCells.clear();
}

farsa::RandomGenerator* CartPoleExperiment::getRng()
{
	return &m_rng;
}

void CartPoleExperiment::setSeed(int s)
{
	m_rng.setSeed(s);
}

}

