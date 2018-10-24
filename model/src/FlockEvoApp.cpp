#include <cmath>
#include <array>
#include <fstream>
#include <iostream>
#include <thread>
#include <ppl.h>
#include "AppSettings.h"
#include "Predator.h"
#include "Prey.h"
#include "timer.hpp"

class FlockEvoApp 
{
public:
  FlockEvoApp();
  void run();

private:
  void simulateGeneration(Predator &predator, std::vector<Prey> tempPrey);
  void createNewGeneration();

  std::vector<Predator> predatorAnimats;

  int generationCount;
  int stepCount;
  int currentPredator;
  int totalHuntCount;

  std::ofstream CSVFile;
};

FlockEvoApp::FlockEvoApp()
{
  rnd_seed((unsigned)time(NULL));
  srand(time(NULL));
}

void FlockEvoApp::run()
{
  // for all iterations
  for (int iterationCount = 0; iterationCount < AppSettings::noOfIterations; iterationCount++)
  {
    // init
    generationCount = 0;
    stepCount = 0;
    currentPredator = 0;
		totalHuntCount = 0;

    // log file
    char fileName[64];
    sprintf_s(fileName, 64, "log%d", iterationCount + 1);
    char extension[8] = ".csv";
    strcat_s(fileName, 64, extension);
    CSVFile = std::ofstream(fileName);

    // create predators
    predatorAnimats.clear();
    for(int i = 0; i < AppSettings::noOfPredatorAnimats; i++)
    {
      predatorAnimats.push_back(Predator(i + 1));
    }

    for (generationCount = 0; generationCount < AppSettings::noOfGenerations; generationCount++)
    {
      if (generationCount != 0)
        createNewGeneration();

      timer tg;

			for (int flockCount = 0; flockCount < AppSettings::noOfFlocks; flockCount++)
			{
				// run noOfSteps for all predators
				std::vector<Prey> tempAnimats;
				for (int j = 0; j < AppSettings::noOfPreyAnimats; j++)
				{
					tempAnimats.push_back(Prey(j));
				}
				concurrency::parallel_for_each(predatorAnimats.begin(), predatorAnimats.end(), [this, &tempAnimats](Predator& pred){
					simulateGeneration(pred, tempAnimats);
				});
			}
	  
			int reportHuntCount = 0;
      for (Predator const& predator : predatorAnimats)
      {
        totalHuntCount += predator.huntCount;
         reportHuntCount = (predator.huntCount - AppSettings::initialHuntCount) / AppSettings::huntFactor;
				CSVFile << iterationCount + 1 << ';';
				CSVFile << generationCount + 1 << ';';
				CSVFile << predator.id << ';';
				CSVFile << (predator.lockOnDistance / (AppSettings::preySize * 2.0f)) << ';';
				CSVFile << (predator.lockOnRadius / (AppSettings::preySize * 2.0f)) << ';';
				CSVFile << predator.distanceForAcceleration << ';';
				CSVFile << predator.velocityMultiplier << ';';
				CSVFile << predator.attackPeriod << ';';
				CSVFile << reportHuntCount;
				CSVFile << std::endl;
			}

      std::cout << iterationCount + 1 << ' ' << generationCount + 1 << ' ';
	  std::cout << tg.elapsed() << 's' << ' '<< reportHuntCount << std::endl;
    }
  }
}

// search pairwise neighbours
void searchNeighbours(std::vector<Prey>& prey)
{
  const float maxSize = std::max(AppSettings::cohesionSize, AppSettings::peripheralitySize);
  const float maxSize2 = maxSize * maxSize;
  auto end = prey.end();
  for (auto a = prey.begin(); a != end; ++a)
  {
    if (!a->isDead)
    {
      for (auto b = a + 1; b != end; ++b)
      {
        if (!b->isDead)
        {
          const glm::vec2 ofs = b->position - a->position;
          const float dist2 = glm::length2(ofs);
          if ((dist2 > 0.0f) && (dist2 < maxSize2))
          {
						float dist = std::sqrt(dist2) - (2 * AppSettings::preySize);
						if (dist < .0f)
							dist = .0f;
						if (dist)
						{
							a->neighbours.emplace_back(ofs, b->getVelocity(), dist);
							b->neighbours.emplace_back(-ofs, a->getVelocity(), dist);
						}
          }
        }
      }
    }
  }
}

void FlockEvoApp::simulateGeneration(Predator &predator, std::vector<Prey> tempPrey)
{
  for (int i = 0; i < AppSettings::noOfSteps; i++)
  {
    // calculate
    // prey
    searchNeighbours(tempPrey);
    for(std::vector<Prey>::iterator p = tempPrey.begin(); p != tempPrey.end(); ++p)
    {
      if (!p->isDead)
        p->calculate(predator);
    }
    // predator
    predator.calculate(tempPrey);

    // update
    // prey
    for(std::vector<Prey>::iterator p = tempPrey.begin(); p != tempPrey.end(); ++p)
    {
      if (!p->isDead)
        p->update(predator, tempPrey);
    }
    // predator
    predator.update();
  }

	predator.reset();
}

void FlockEvoApp::createNewGeneration()
{
  //genetic stuff
  std::vector<Predator> oldPredators = predatorAnimats;
  predatorAnimats.clear();
  for(int i = 0; i < AppSettings::noOfPredatorAnimats; i++)
  {
    // first select 2 random parents
    int randomIndex = randomInt(0, AppSettings::noOfPredatorAnimats - 1);
    Predator firstParent = oldPredators.at(randomIndex);
    randomIndex = randomInt(0, AppSettings::noOfPredatorAnimats - 1);
    Predator secondParent = oldPredators.at(randomIndex);

    // get parents
    if (totalHuntCount > 1)
    {
      int currentHuntCount = 0;
      int randomHuntCount = randomInt(1, totalHuntCount);
      // first parent
      for (int j = 0; j < AppSettings::noOfPredatorAnimats; j++)
      {
        firstParent = oldPredators.at(j);
        currentHuntCount += firstParent.huntCount;
        if (currentHuntCount >= randomHuntCount)
          break;
      }

      currentHuntCount = 0;
      randomHuntCount = randomInt(1, totalHuntCount);
      // second parent
      for (int j = 0; j < AppSettings::noOfPredatorAnimats; j++)
      {
        secondParent = oldPredators.at(j);
        currentHuntCount += secondParent.huntCount;
        if (currentHuntCount >= randomHuntCount)
          break;
      }
    }

	int aP;
	float vM;
	float dA;

	float random = randomFloat(.0f, 1.0f);
	if (random < 0.5f) aP = firstParent.attackPeriod;
	else aP = secondParent.attackPeriod;

	random = randomFloat(.0f, 1.0f);
	if (random < 0.5f) vM = firstParent.velocityMultiplier;
	else vM = secondParent.velocityMultiplier;

	random = randomFloat(.0f, 1.0f);
	if (random < 0.5f) dA = firstParent.distanceForAcceleration;
	else dA = secondParent.distanceForAcceleration;

	// mutation

	random = randomFloat(.0f, 1.0f);
	if (random < AppSettings::mutationRate){
		float flip = randomFloat(.0f, 1.0f);
		if (flip > 0.5f) aP = aP + (aP * AppSettings::mutationFactor);
		else aP = aP - (aP * AppSettings::mutationFactor);
	}
	if (aP < 0) aP = 0;

	random = randomFloat(.0f, 1.0f);
	if (random < AppSettings::mutationRate){
		float flip = randomFloat(.0f, 1.0f);
		if (flip > 0.5f) vM += AppSettings::mutationFactor;
		else vM -= AppSettings::mutationFactor;
	}
	if (vM < 1.0f) vM = 1.0f;

	random = randomFloat(.0f, 1.0f);
	if (random < AppSettings::mutationRate){
		float flip = randomFloat(.0f, 1.0f);
		if (flip > 0.5f) dA = dA + (dA * AppSettings::mutationFactor);
		else dA = dA - (dA * AppSettings::mutationFactor);;
	}
	if (dA < 0) dA = 0;

    // crossover - lock-on distance from first parent, lock-on radius from second
    float parentLockOnDistance = firstParent.lockOnDistance;
    float parentLockOnRadius = secondParent.lockOnRadius;

    // mutation
    random = randomFloat(.0f, 1.0f);
    if (random < AppSettings::mutationRate)
    {
      float flip = randomFloat(.0f, 1.0f);
      if (flip > 0.5f)
        parentLockOnDistance = parentLockOnDistance + (parentLockOnDistance * AppSettings::mutationFactor);
      else
        parentLockOnDistance = parentLockOnDistance - (parentLockOnDistance * AppSettings::mutationFactor);

      if (parentLockOnDistance < 0)
        parentLockOnDistance = 1.0f;
      else if (parentLockOnDistance > AppSettings::huntSize)
        parentLockOnDistance = AppSettings::huntSize;
    }

      random = randomFloat(.0f, 1.0f);
      if (random < AppSettings::mutationRate)
      {
        float flip = randomFloat(.0f, 1.0f);
        if (flip > 0.5f)
          parentLockOnRadius = parentLockOnRadius + (parentLockOnRadius * AppSettings::mutationFactor);
        else
          parentLockOnRadius = parentLockOnRadius - (parentLockOnRadius * AppSettings::mutationFactor);

        if (parentLockOnRadius < 0)
          parentLockOnRadius = 1.0f;
        else if (parentLockOnRadius > AppSettings::huntSize)
          parentLockOnRadius = AppSettings::huntSize;
      }

    predatorAnimats.push_back(Predator(i + 1, parentLockOnDistance, parentLockOnRadius, dA, vM, aP));
  }

  currentPredator = 0;
  totalHuntCount = 0;
}

int main()
{
  FlockEvoApp app;    
  app.run();
}