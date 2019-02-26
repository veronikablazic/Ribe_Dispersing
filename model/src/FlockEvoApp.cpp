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
  void simulateGeneration(Predator &predator, std::vector<Prey> tempPrey, int index);
  void createNewGeneration();

  std::vector<Predator> predatorAnimats;
  std::vector<std::vector<Prey>> preyAnimats;

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

	preyAnimats.clear();

	if (PREY_EVOLUTION == 1) {
		for (int j = 0; j < AppSettings::noOfPreyAnimats; j++) {
			std::vector<Prey> tempPrey;
			for (int i = 0; i < AppSettings::noOfPreyAnimats; i++) {
				tempPrey.push_back(Prey(i));
			}
			preyAnimats.push_back(tempPrey);
		}
	}

    for (generationCount = 0; generationCount < AppSettings::noOfGenerations; generationCount++)
    {
      if (generationCount != 0)
        createNewGeneration();

      timer tg;

		#if(PREY_EVOLUTION == 0)
			  for (int flockCount = 0; flockCount < AppSettings::noOfFlocks; flockCount++) {

				  std::vector<Prey> tempPrey;
				  for (int j = 0; j < AppSettings::noOfPreyAnimats; j++) {
					  tempPrey.push_back(Prey(j));
				  }

				  concurrency::parallel_for_each(predatorAnimats.begin(), predatorAnimats.end(), [this, &tempPrey](Predator& pred){
					  simulateGeneration(pred, tempPrey, 0);
				  });

			  }
		#endif
			  // part 2
		#if(PREY_EVOLUTION == 1)

			  // predator nas tu ne zanima
			  for (int flockCount = 0; flockCount < AppSettings::noOfFlocks; flockCount++) {

				  int n = preyAnimats.size();

				  concurrency::parallel_for(0, n, [&](int index) {
					  Predator tempPred = Predator(0);
					  simulateGeneration(tempPred, preyAnimats[index], index);
				  });
			  }

		#endif
	  
			char buffer[128];
			int reportHuntCount = 0;
#if (PREY_EVOLUTION == 1)
			int count = 0;
			for (std::vector<Prey> prey : preyAnimats)
			{
				count++;
				deadCount = 0;
				for (int k = 0; k < 100; k++) {
					if (prey[k].isDead) deadCount++;
					else totalEnergy += prey[k].energy;
				}

				// povprecna preostala energija
				totalEnergy = totalEnergy / (prey.size() - deadCount);

				CSVFile << iterationCount + 1 << ';';
				CSVFile << generationCount + 1 << ';';
				CSVFile << count << ';';
				CSVFile << prey[0].wb << ';';
				CSVFile << prey[0].wr << ';';
				CSVFile << prey[0].we << ';';
				CSVFile << prey[0].cn0 << ';';
				CSVFile << prey[0].cn1 << ';';
				CSVFile << prey[0].cn2 << ';';
				CSVFile << prey[0].cn3 << ';';
				CSVFile << prey[0].cr0 << ';';
				CSVFile << prey[0].cr1 << ';';
				CSVFile << prey[0].selfishTime << ';';
				CSVFile << prey[0].selfishTime2 << ';';
				CSVFile << prey[0].selfishEscapeDistance << ';';
				CSVFile << prey[0].selfishProbability << ';';
				CSVFile << prey[0].dot_pow << ';';
				CSVFile << prey[0].dotTreshold << ';';
				CSVFile << prey[0].v_b << ';';
				CSVFile << prey[0].v_r << ';';
				CSVFile << totalEnergy << ';';
				CSVFile << deadCount;
				CSVFile << std::endl;
			}

			sprintf_s(buffer, "%d %d %d %gs\n", iterationCount + 1, generationCount + 1, deadCount, tg.elapsed());
			std::cout << buffer;
#endif

#if (PREY_EVOLUTION == 0)
			reportHuntCount = 0;
			for (Predator const& predator : predatorAnimats)
			{
				totalHuntCount += predator.huntCount;

				// part 1
				if (PRED_ENERGY_PARAMS == 0) {
					//int reportHuntCount = (predator.huntCount - AppSettings::initialHuntCount) / AppSettings::huntFactor;
					//reportHuntCount += predator.huntCount;
					reportHuntCount = predator.huntCount;
					CSVFile << iterationCount + 1 << ';';
					CSVFile << generationCount + 1 << ';';
					CSVFile << predator.id << ';';
					CSVFile << (predator.lockOnDistance / (AppSettings::preySize * 2.0f)) << ';';
					CSVFile << (predator.lockOnRadius / (AppSettings::preySize * 2.0f)) << ';';
					CSVFile << predator.huntCount;
					CSVFile << std::endl;
				}
				// part 3
				else {
					//int reportHuntCount = (predator.huntCount - AppSettings::initialHuntCount) / AppSettings::huntFactor;
					//reportHuntCount += predator.huntCount;
					reportHuntCount = predator.huntCount;
					CSVFile << iterationCount + 1 << ';';
					CSVFile << generationCount + 1 << ';';
					CSVFile << predator.id << ';';
					CSVFile << (predator.lockOnDistance / (AppSettings::preySize * 2.0f)) << ';';
					CSVFile << (predator.lockOnRadius / (AppSettings::preySize * 2.0f)) << ';';
					CSVFile << predator.distanceForAcceleration << ';';
					CSVFile << predator.wanderingTime << ';';
					CSVFile << predator.attackPeriod << ';';
					CSVFile << predator.wb << ';';
					CSVFile << predator.wr << ';';
					CSVFile << predator.cn0 << ';';
					CSVFile << predator.cn1 << ';';
					CSVFile << predator.cn2 << ';';
					CSVFile << predator.cn3 << ';';
					CSVFile << predator.cr0 << ';';
					CSVFile << predator.cr1 << ';';
					CSVFile << predator.v_b << ';';
					CSVFile << predator.v_r << ';';
					CSVFile << predator.energy << ';';
					CSVFile << predator.huntCount;
					CSVFile << std::endl;
				}

				//sprintf_s(buffer, "%d;%d;%d;%.2f;%.2f;%.2f;%d\n", iterationCount + 1, generationCount + 1, predator.id, predator.nearestWeight, predator.peripheralWeight, predator.centreWeight, reportHuntCount);
				//CSVFile << buffer;
			}

			sprintf_s(buffer, "%d %d %d %g s\n", iterationCount + 1, generationCount + 1, reportHuntCount, tg.elapsed());
			std::cout << buffer;
#endif
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

void FlockEvoApp::simulateGeneration(Predator &predator, std::vector<Prey> tempPrey, int index)
{
	for (int i = 0; i < AppSettings::noOfSteps; i++) {

		// calculate
		// prey
		searchNeighbours(tempPrey);
		for (std::vector<Prey>::iterator p = tempPrey.begin(); p != tempPrey.end(); ++p)
		{
			if (!p->isDead)
				p->calculate(predator);
		}
		// predator
		predator.calculate(tempPrey);

		// update
		// prey
		for (std::vector<Prey>::iterator p = tempPrey.begin(); p != tempPrey.end(); ++p)
		{
			if (!p->isDead)
				p->update(predator, tempPrey);
		}
		// predator
		predator.update(tempPrey);

		if (predator.isExhausted) break;
	}

	if (PREY_EVOLUTION == 1) {
		preyAnimats[index] = tempPrey;
		for (int j = 0; j < tempPrey.size(); j++) preyAnimats[index][j].reset();
	}
	predator.reset();
}

void FlockEvoApp::createNewGeneration()
{
#if(PREY_EVOLUTION == 0)
	//genetic stuff
	std::vector<Predator> oldPredators = predatorAnimats;
	predatorAnimats.clear();
	for (int i = 0; i < AppSettings::noOfPredatorAnimats; i++)
	{
		// first select 2 random parents
		int randomIndex = randomInt(0, AppSettings::noOfPredatorAnimats - 1);
		Predator firstParent = oldPredators.at(randomIndex);
		randomIndex = randomInt(0, AppSettings::noOfPredatorAnimats - 1);
		Predator secondParent = oldPredators.at(randomIndex);

		// get parents
		if (PRED_ENERGY_PARAMS == 0) {
			if (totalHuntCount > 1) {
				int currentHuntCount = 0;
				int randomHuntCount = randomInt(1, totalHuntCount);
				// first parent
				for (int j = 0; j < AppSettings::noOfPredatorAnimats; j++) {
					firstParent = oldPredators.at(j);
					currentHuntCount += firstParent.huntCount;
					if (currentHuntCount >= randomHuntCount) break;
				}

				currentHuntCount = 0;
				randomHuntCount = randomInt(1, totalHuntCount);
				// second parent
				for (int j = 0; j < AppSettings::noOfPredatorAnimats; j++){
					secondParent = oldPredators.at(j);
					currentHuntCount += secondParent.huntCount;
					if (currentHuntCount >= randomHuntCount)
						break;
				}
			}
		}
		else {

			totalHuntCount = 20;

			if (totalHuntCount > 0) {

				float maxDead = 0;
				float maxEnergy = 0;
				for (int j = 0; j < oldPredators.size(); j++) {
					if (oldPredators[j].huntCount > maxDead) maxDead = oldPredators[j].huntCount;
					if (oldPredators[j].energy > maxEnergy) maxEnergy = oldPredators[j].energy;
				}

				float firstFitness = 0.0f;
				int firstParentInt = 0;
				float secondFitness = 0.0f;
				int secondParentInt = 0;

				// vzamemo 2 najboljsa
				for (int j = 0; j < oldPredators.size(); j++) {
					float fitness = (float)(oldPredators[j].huntCount / maxDead) + oldPredators[j].step / 600.0f + oldPredators[j].energy / maxEnergy;
					if (fitness > firstFitness) {
						secondFitness = firstFitness;
						secondParentInt = firstParentInt;
						firstFitness = fitness;
						firstParentInt = j;
					}
					else if (fitness > secondFitness) {
						secondFitness = fitness;
						secondParentInt = j;
					}
				}

				firstParent = oldPredators.at(firstParentInt);
				secondParent = oldPredators.at(secondParentInt);
			}
		}

		float random;

		float wb1 = 0;
		float wr1 = 0;
		float we1 = 0;
		float cn01 = 0;
		float cn11 = 0;
		float cn21 = 0;
		float cn31 = 0;
		float cr01 = 0;
		float cr11 = 0;
		float v_b1 = 0;
		float v_r1 = 0;

		int aP = 0;
		float dA = 0;
		int wT = 0;

		if (PRED_ENERGY_PARAMS == 1) {

			// coinflip crossover

			float random = randomFloat(0.0f, 1.0f);
			if (random < 0.5f) wb1 = firstParent.wb;
			else wb1 = secondParent.wb;

			random = randomFloat(0.0f, 1.0f);
			if (random < 0.5f) wr1 = firstParent.wr;
			else wr1 = secondParent.wr;

			random = randomFloat(0.0f, 1.0f);
			if (random < 0.5f) we1 = firstParent.we;
			else we1 = secondParent.we;

			random = randomFloat(0.0f, 1.0f);
			if (random < 0.5f) cn01 = firstParent.cn0;
			else cn01 = secondParent.cn0;

			random = randomFloat(0.0f, 1.0f);
			if (random < 0.5f) cn11 = firstParent.cn1;
			else cn11 = secondParent.cn1;

			random = randomFloat(0.0f, 1.0f);
			if (random < 0.5f) cn21 = firstParent.cn2;
			else cn21 = secondParent.cn2;

			random = randomFloat(0.0f, 1.0f);
			if (random < 0.5f) cn31 = firstParent.cn3;
			else cn31 = secondParent.cn3;

			random = randomFloat(0.0f, 1.0f);
			if (random < 0.5f) cr01 = firstParent.cr0;
			else cr01 = secondParent.cr0;

			random = randomFloat(0.0f, 1.0f);
			if (random < 0.5f) cr11 = firstParent.cr1;
			else cr11 = secondParent.cr1;

			random = randomFloat(0.0f, 1.0f);
			if (random < 0.5f) wT = firstParent.wanderingTime;
			else wT = secondParent.wanderingTime;

			random = randomFloat(.0f, 1.0f);
			if (random < 0.5f) aP = firstParent.attackPeriod;
			else aP = secondParent.attackPeriod;

			random = randomFloat(.0f, 1.0f);
			if (random < 0.5f) dA = firstParent.distanceForAcceleration;
			else dA = secondParent.distanceForAcceleration;

			random = randomFloat(0.0f, 1.0f);
			if (random < 0.5f) v_r1 = firstParent.v_r;
			else v_r1 = secondParent.v_r;

			random = randomFloat(0.0f, 1.0f);
			if (random < 0.5f) v_b1 = firstParent.v_b;
			else v_b1 = secondParent.v_b;

			// mutation
			random = randomFloat(0.0f, 1.0f);
			if (random < AppSettings::mutationRate){
				float flip = randomFloat(.0f, 1.0f);
				if (flip > 0.5f) wb1 += AppSettings::mutationFactor;
				else wb1 -= AppSettings::mutationFactor;
			}
			wb1 = glm::clamp(wb1, 0.0f, 1.0f);

			random = randomFloat(0.0f, 1.0f);
			if (random < AppSettings::mutationRate){
				float flip = randomFloat(.0f, 1.0f);
				if (flip > 0.5f) wr1 += AppSettings::mutationFactor;
				else wr1 -= AppSettings::mutationFactor;
			}
			wr1 = glm::clamp(wr1, 0.0f, 1.0f);

			random = randomFloat(0.0f, 1.0f);
			if (random < AppSettings::mutationRate){
				float flip = randomFloat(.0f, 1.0f);
				if (flip > 0.5f) we1 += AppSettings::mutationFactor;
				else we1 -= AppSettings::mutationFactor;
			}
			we1 = glm::clamp(we1, 0.0f, 1.0f);

			random = randomFloat(0.0f, 1.0f);
			if (random < AppSettings::mutationRate){
				float flip = randomFloat(.0f, 1.0f);
				if (flip > 0.5f) cn01 += AppSettings::mutationFactor;
				else cn01 -= AppSettings::mutationFactor;
			}
			cn01 = glm::clamp(cn01, AppSettings::minPredatorVelocity, AppSettings::maxPredatorVelocity);

			random = randomFloat(0.0f, 1.0f);
			if (random < AppSettings::mutationRate){
				float flip = randomFloat(.0f, 1.0f);
				if (flip > 0.5f) cn11 += AppSettings::mutationFactor;
				else cn11 -= AppSettings::mutationFactor;
			}
			cn11 = glm::clamp(cn11, 0.0f, 50.0f);

			random = randomFloat(0.0f, 1.0f);
			if (random < AppSettings::mutationRate){
				float flip = randomFloat(.0f, 1.0f);
				if (flip > 0.5f) cn21 += AppSettings::mutationFactor;
				else cn21 -= AppSettings::mutationFactor;
			}
			cn21 = glm::clamp(cn21, AppSettings::minPredatorVelocity, AppSettings::maxPredatorVelocity);

			random = randomFloat(0.0f, 1.0f);
			if (random < AppSettings::mutationRate){
				float flip = randomFloat(.0f, 1.0f);
				if (flip > 0.5f) cn31 += AppSettings::mutationFactor;
				else cn31 -= AppSettings::mutationFactor;
			}
			cn31 = glm::clamp(cn31, 0.0f, 50.0f);

			random = randomFloat(0.0f, 1.0f);
			if (random < AppSettings::mutationRate){
				float flip = randomFloat(.0f, 1.0f);
				if (flip > 0.5f) cr01 += AppSettings::mutationFactor;
				else cr01 -= AppSettings::mutationFactor;
			}
			cr01 = glm::clamp(cr01, 0.0f, 1.0f);

			random = randomFloat(0.0f, 1.0f);
			if (random < AppSettings::mutationRate){
				float flip = randomFloat(.0f, 1.0f);
				if (flip > 0.5f) cr11 += AppSettings::mutationFactor;
				else cr11 -= AppSettings::mutationFactor;
			}
			cr11 = glm::clamp(cr11, 0.0f, 50.0f);

			random = randomFloat(0.0f, 1.0f);
			if (random < AppSettings::mutationRate){
				float flip = randomFloat(.0f, 1.0f);
				if (flip > 0.5f) v_r1 += AppSettings::mutationFactor;
				else v_r1 -= AppSettings::mutationFactor;
			}
			v_r1 = glm::clamp(v_r1, AppSettings::minPredatorVelocity, AppSettings::cruisingSpeed);

			random = randomFloat(0.0f, 1.0f);
			if (random < AppSettings::mutationRate){
				float flip = randomFloat(.0f, 1.0f);
				if (flip > 0.5f) v_b1 += AppSettings::mutationFactor;
				else v_b1 -= AppSettings::mutationFactor;
			}
			v_b1 = glm::clamp(v_b1, AppSettings::cruisingSpeed, AppSettings::maxPredatorVelocity);

			random = randomFloat(.0f, 1.0f);
			if (random < AppSettings::mutationRate){
				float flip = randomFloat(.0f, 1.0f);
				if (flip > 0.5f) aP += AppSettings::mutationFactorInt;
				else aP -= AppSettings::mutationFactorInt;
			}
			if (aP < 0) aP = 0;
			if (aP > AppSettings::noOfSteps) aP = AppSettings::noOfSteps;

			random = randomFloat(.0f, 1.0f);
			if (random < AppSettings::mutationRate){
				float flip = randomFloat(.0f, 1.0f);
				if (flip > 0.5f) dA += AppSettings::mutationFactor;
				else dA -= AppSettings::mutationFactor;
			}
			if (dA < 0) dA = 0;
			if (dA > AppSettings::huntSize) dA = AppSettings::huntSize;

			random = randomFloat(.0f, 1.0f);
			if (random < AppSettings::mutationRate){
				float flip = randomFloat(.0f, 1.0f);
				if (flip > 0.5f) wT += AppSettings::mutationFactorInt;
				else wT -= AppSettings::mutationFactorInt;
			}
			if (wT < 0) wT = 0;
			if (wT > AppSettings::noOfSteps) wT = AppSettings::noOfSteps;
		}

		// coinflip crossover
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
		predatorAnimats.push_back(Predator(i + 1, parentLockOnDistance, parentLockOnRadius, dA, aP, wT, wb1, wr1, we1, cn01, cn11, cn21, cn31, cr01, cr11, v_b1, v_r1));
	}

	currentPredator = 0;
	totalHuntCount = 0;
#endif

#if (PREY_EVOLUTION == 1)

	std::vector<std::vector<Prey>> oldPrey = preyAnimats;
	preyAnimats.clear();

	float totalPreyFitness = 0.0f;
	std::vector<int> noAlive;
	std::vector<float> avgEnergy;
	float maxEnergy = 0;
	std::vector<float> preyFitness;

	for (int i = 0; i < oldPrey.size(); i++) {
		float aliveCount = 0;
		float totalEnergy = 0;
		for (int j = 0; j < AppSettings::noOfPreyAnimats; j++) {
			if (!oldPrey[i][j].isDead) {
				aliveCount++;
				totalEnergy += oldPrey[i][j].energy;
			}
		}
		noAlive.push_back(aliveCount);
		avgEnergy.push_back(totalEnergy / aliveCount);
		if (totalEnergy / aliveCount > maxEnergy) maxEnergy = totalEnergy / aliveCount;
	}

	// get fitness for each flock
	for (int i = 0; i < oldPrey.size(); i++) {
		// energija pomnozena z utezjo, ker je manj pomembna od ulova
		preyFitness.push_back((noAlive[i] / 100.0f) + (avgEnergy[i] / maxEnergy));
		totalPreyFitness += preyFitness[i];
	}

	for (int i = 0; i < 100; i++) {

		int randomIndex = randomInt(0, oldPrey.size() - 1);
		std::vector<Prey> firstParent = oldPrey.at(randomIndex);

		randomIndex = randomInt(0, oldPrey.size() - 1);
		std::vector<Prey> secondParent = oldPrey.at(randomIndex);

		if (totalPreyFitness > 0) {

			float firstFitness = 0.0f;
			int firstParentInt = 0;
			float secondFitness = 0.0f;
			int secondParentInt = 0;

			// vzamemo 2 najboljsa
			for (int j = 0; j < preyFitness.size(); j++) {
				if (preyFitness[j] > firstFitness) {
					secondFitness = firstFitness;
					secondParentInt = firstParentInt;
					firstFitness = preyFitness[j];
					firstParentInt = j;
				}
				else if (preyFitness[j] > secondFitness) {
					secondFitness = preyFitness[j];
					secondParentInt = j;
				}
			}

			firstParent = oldPrey.at(firstParentInt);
			secondParent = oldPrey.at(secondParentInt);

			/*float currentFitness = 0;
			float randomFitness = randomFloat(0.f, totalPreyFitness);

			for (int j = 0; j < oldPrey.size(); j++) {
			firstParent = oldPrey.at(j);
			currentFitness += preyFitness[j];
			if (currentFitness > randomFitness) break;
			}

			currentFitness = 0;
			randomFitness = randomFloat(0.f, totalPreyFitness);

			for (int j = 0; j < oldPrey.size(); j++) {
			secondParent = oldPrey.at(j);
			currentFitness += preyFitness[j];
			if (currentFitness >= randomFitness) break;
			}*/
		}

		float wb1, wr1, we1;
		float cn01, cn11, cn21, cn31;
		float cr01, cr11;
		int selfishTime1, selfishTime21;
		float selfishEscapeDistance1, selfishProbability1;
		float dot_pow1, dotTreshold1;
		float v_b1, v_r1;

		// coinflip crossover

		float random = randomFloat(0.0f, 1.0f);
		if (random < 0.5f) wb1 = firstParent[0].wb;
		else wb1 = secondParent[0].wb;

		random = randomFloat(0.0f, 1.0f);
		if (random < 0.5f) wr1 = firstParent[0].wr;
		else wr1 = secondParent[0].wr;

		random = randomFloat(0.0f, 1.0f);
		if (random < 0.5f) we1 = firstParent[0].we;
		else we1 = secondParent[0].we;

		random = randomFloat(0.0f, 1.0f);
		if (random < 0.5f) cn01 = firstParent[0].cn0;
		else cn01 = secondParent[0].cn0;

		random = randomFloat(0.0f, 1.0f);
		if (random < 0.5f) cn11 = firstParent[0].cn1;
		else cn11 = secondParent[0].cn1;

		random = randomFloat(0.0f, 1.0f);
		if (random < 0.5f) cn21 = firstParent[0].cn2;
		else cn21 = secondParent[0].cn2;

		random = randomFloat(0.0f, 1.0f);
		if (random < 0.5f) cn31 = firstParent[0].cn3;
		else cn31 = secondParent[0].cn3;

		random = randomFloat(0.0f, 1.0f);
		if (random < 0.5f) cr01 = firstParent[0].cr0;
		else cr01 = secondParent[0].cr0;

		random = randomFloat(0.0f, 1.0f);
		if (random < 0.5f) cr11 = firstParent[0].cr1;
		else cr11 = secondParent[0].cr1;

		random = randomFloat(0.0f, 1.0f);
		if (random < 0.5f) selfishTime1 = firstParent[0].selfishTime;
		else selfishTime1 = secondParent[0].selfishTime;

		random = randomFloat(0.0f, 1.0f);
		if (random < 0.5f) selfishTime21 = firstParent[0].selfishTime2;
		else selfishTime21 = secondParent[0].selfishTime2;

		random = randomFloat(0.0f, 1.0f);
		if (random < 0.5f) selfishEscapeDistance1 = firstParent[0].selfishEscapeDistance;
		else selfishEscapeDistance1 = secondParent[0].selfishEscapeDistance;

		random = randomFloat(0.0f, 1.0f);
		if (random < 0.5f) selfishProbability1 = firstParent[0].selfishProbability;
		else selfishProbability1 = secondParent[0].selfishProbability;

		random = randomFloat(0.0f, 1.0f);
		if (random < 0.5f) dot_pow1 = firstParent[0].dot_pow;
		else dot_pow1 = secondParent[0].dot_pow;

		random = randomFloat(0.0f, 1.0f);
		if (random < 0.5f) dotTreshold1 = firstParent[0].dotTreshold;
		else dotTreshold1 = secondParent[0].dotTreshold;

		random = randomFloat(0.0f, 1.0f);
		if (random < 0.5f) v_r1 = firstParent[0].v_r;
		else v_r1 = secondParent[0].v_r;

		random = randomFloat(0.0f, 1.0f);
		if (random < 0.5f) v_b1 = firstParent[0].v_b;
		else v_b1 = secondParent[0].v_b;

		// mutation
		random = randomFloat(0.0f, 1.0f);
		if (random < AppSettings::mutationRate){
			float flip = randomFloat(.0f, 1.0f);
			if (flip > 0.5f) wb1 += AppSettings::mutationFactor;
			else wb1 -= AppSettings::mutationFactor;
		}
		wb1 = glm::clamp(wb1, 0.0f, 1.0f);

		random = randomFloat(0.0f, 1.0f);
		if (random < AppSettings::mutationRate){
			float flip = randomFloat(.0f, 1.0f);
			if (flip > 0.5f) wr1 += AppSettings::mutationFactor;
			else wr1 -= AppSettings::mutationFactor;
		}
		wr1 = glm::clamp(wr1, 0.0f, 1.0f);

		random = randomFloat(0.0f, 1.0f);
		if (random < AppSettings::mutationRate){
			float flip = randomFloat(.0f, 1.0f);
			if (flip > 0.5f) we1 += AppSettings::mutationFactor;
			else we1 -= AppSettings::mutationFactor;
		}
		we1 = glm::clamp(we1, 0.0f, 1.0f);

		random = randomFloat(0.0f, 1.0f);
		if (random < AppSettings::mutationRate){
			float flip = randomFloat(.0f, 1.0f);
			if (flip > 0.5f) cn01 += AppSettings::mutationFactor;
			else cn01 -= AppSettings::mutationFactor;
		}
		cn01 = glm::clamp(cn01, AppSettings::minPreyVelocity, AppSettings::maxPreyVelocity);

		random = randomFloat(0.0f, 1.0f);
		if (random < AppSettings::mutationRate){
			float flip = randomFloat(.0f, 1.0f);
			if (flip > 0.5f) cn11 += AppSettings::mutationFactor;
			else cn11 -= AppSettings::mutationFactor;
		}
		cn11 = glm::clamp(cn11, 0.0f, 50.0f);

		random = randomFloat(0.0f, 1.0f);
		if (random < AppSettings::mutationRate){
			float flip = randomFloat(.0f, 1.0f);
			if (flip > 0.5f) cn21 += AppSettings::mutationFactor;
			else cn21 -= AppSettings::mutationFactor;
		}
		cn21 = glm::clamp(cn21, AppSettings::minPreyVelocity, AppSettings::maxPreyVelocity);

		random = randomFloat(0.0f, 1.0f);
		if (random < AppSettings::mutationRate){
			float flip = randomFloat(.0f, 1.0f);
			if (flip > 0.5f) cn31 += AppSettings::mutationFactor;
			else cn31 -= AppSettings::mutationFactor;
		}
		cn31 = glm::clamp(cn31, 0.0f, 50.0f);

		random = randomFloat(0.0f, 1.0f);
		if (random < AppSettings::mutationRate){
			float flip = randomFloat(.0f, 1.0f);
			if (flip > 0.5f) cr01 += AppSettings::mutationFactor;
			else cr01 -= AppSettings::mutationFactor;
		}
		cr01 = glm::clamp(cr01, 0.0f, 1.0f);

		random = randomFloat(0.0f, 1.0f);
		if (random < AppSettings::mutationRate){
			float flip = randomFloat(.0f, 1.0f);
			if (flip > 0.5f) cr11 += AppSettings::mutationFactor;
			else cr11 -= AppSettings::mutationFactor;
		}
		cr11 = glm::clamp(cr11, 0.0f, 50.0f);

		random = randomFloat(0.0f, 1.0f);
		if (random < AppSettings::mutationRate){
			float flip = randomFloat(.0f, 1.0f);
			if (flip > 0.5f) selfishTime1 += AppSettings::mutationFactor;
			else selfishTime1 -= AppSettings::mutationFactor;
		}
		selfishTime1 = glm::clamp(selfishTime1, 1, 600);

		random = randomFloat(0.0f, 1.0f);
		if (random < AppSettings::mutationRate){
			float flip = randomFloat(.0f, 1.0f);
			if (flip > 0.5f) selfishTime21 += AppSettings::mutationFactor;
			else selfishTime21 -= AppSettings::mutationFactor;
		}
		selfishTime21 = glm::clamp(selfishTime21, 0, 600);

		random = randomFloat(0.0f, 1.0f);
		if (random < AppSettings::mutationRate){
			float flip = randomFloat(.0f, 1.0f);
			if (flip > 0.5f) selfishEscapeDistance1 += AppSettings::mutationFactor;
			else selfishEscapeDistance1 -= AppSettings::mutationFactor;
		}
		selfishEscapeDistance1 = glm::clamp(selfishEscapeDistance1, 0.0f, AppSettings::huntSize);

		random = randomFloat(0.0f, 1.0f);
		if (random < AppSettings::mutationRate){
			float flip = randomFloat(.0f, 1.0f);
			if (flip > 0.5f) selfishProbability1 += AppSettings::mutationFactor;
			else selfishProbability1 -= AppSettings::mutationFactor;
		}
		selfishProbability1 = glm::clamp(selfishProbability1, 0.0f, 1.0f);

		random = randomFloat(0.0f, 1.0f);
		if (random < AppSettings::mutationRate){
			float flip = randomFloat(.0f, 1.0f);
			if (flip > 0.5f) dot_pow1 += AppSettings::mutationFactor;
			else dot_pow1 -= AppSettings::mutationFactor;
		}
		dot_pow1 = glm::clamp(dot_pow1, 0.0f, 50.0f);

		random = randomFloat(0.0f, 1.0f);
		if (random < AppSettings::mutationRate){
			float flip = randomFloat(.0f, 1.0f);
			if (flip > 0.5f) dotTreshold1 += AppSettings::mutationFactor;
			else dotTreshold1 -= AppSettings::mutationFactor;
		}
		dotTreshold1 = glm::clamp(dotTreshold1, 0.0f, 1.0f);

		random = randomFloat(0.0f, 1.0f);
		if (random < AppSettings::mutationRate){
			float flip = randomFloat(.0f, 1.0f);
			if (flip > 0.5f) v_r1 += AppSettings::mutationFactor;
			else v_r1 -= AppSettings::mutationFactor;
		}
		v_r1 = glm::clamp(v_r1, AppSettings::minPreyVelocity, AppSettings::cruisingPrey);

		random = randomFloat(0.0f, 1.0f);
		if (random < AppSettings::mutationRate){
			float flip = randomFloat(.0f, 1.0f);
			if (flip > 0.5f) v_b1 += AppSettings::mutationFactor;
			else v_b1 -= AppSettings::mutationFactor;
		}
		v_b1 = glm::clamp(v_b1, AppSettings::cruisingPrey, AppSettings::maxPreyVelocity);

		std::vector<Prey> tempPrey;
		for (int j = 0; j < AppSettings::noOfPreyAnimats; j++){
			tempPrey.push_back(Prey(j, wb1, wr1, we1, cn01, cn11, cn21, cn31, cr01, cr11, selfishTime1, selfishTime21, selfishEscapeDistance1, selfishProbability1, dot_pow1, dotTreshold1, v_b1, v_r1));
		}
		preyAnimats.push_back(tempPrey);

	}

#endif
}

int main()
{
  FlockEvoApp app;    
  app.run();
}