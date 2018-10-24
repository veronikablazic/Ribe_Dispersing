#include <cassert>
#include "AppSettings.h"
#include "Predator.h"
#include "Prey.h"
#include <algorithm>

Predator::Predator() {}

Predator::Predator(int animatID)
{
  // acceleration
  acceleration = glm::vec2(.0f, .0f);

  // position
  position = glm::vec2(AppSettings::worldCentre.x, AppSettings::worldCentre.y + (AppSettings::preySize * 2.0f * 200.0f));

  // speed and heading
  heading = glm::vec2(.0f, -1.0f);
  speed = randomFloat(AppSettings::minPredatorVelocity, AppSettings::maxPredatorVelocity);

  id = animatID;
  target = -1;
  central = -1;
  huntCount = 0;
  lockOn = false;
  handling = false;
  handlingTimer = 0;

  // random lock-on distances
  lockOnDistance = randomFloat(1.0f, AppSettings::huntSize);
  lockOnRadius = randomFloat(1.0f, AppSettings::huntSize);

  // new 

  energy = 1;
  regenerationGain = 30.0f / 100000.0f; //regeneration in 30'

  attackZoneAngle = (std::_Pi / 6);
  angle_t = 0;

  distFromTarget = 1000000;

  distanceForAcceleration = randomFloat(1.0f, AppSettings::huntSize);
  velocityMultiplier = randomFloat(1.0f, 3.0f);
  attackPeriod = randomFloat(1.0f, AppSettings::noOfSteps);
  currentAttackTime = 0;
}

Predator::Predator(int animatID, float parentLockOnDistance, float parentlockOnRadius, float dA, float vM, int aP)
{
  // acceleration
  acceleration = glm::vec2(.0f, .0f);

  // position
  position = glm::vec2(AppSettings::worldCentre.x, AppSettings::worldCentre.y + (AppSettings::preySize * 2.0f * 200.0f));

  // speed and heading
  heading = glm::vec2(.0f, -1.0f);
  speed = randomFloat(AppSettings::minPredatorVelocity, AppSettings::maxPredatorVelocity);

  id = animatID;
  target = -1;
  central = -1;
  huntCount = AppSettings::initialHuntCount;
  lockOn = false;
  handling = false;
  handlingTimer = 0;

  // lock-on distances from parents
  lockOnDistance = parentLockOnDistance;
  lockOnRadius = parentlockOnRadius;

  energy = 1;
  regenerationGain = 30.0f / 100000.0f; //regeneration in 30'

  attackZoneAngle = (std::_Pi / 6);
  angle_t = 0;

  distFromTarget = 1000000;

  distanceForAcceleration = dA;
  velocityMultiplier = vM;
  attackPeriod = aP;
  currentAttackTime = 0;
}

void Predator::reset()
{
	// acceleration
	acceleration = glm::vec2(.0f, .0f);

	// position
	position = glm::vec2(AppSettings::worldCentre.x, AppSettings::worldCentre.y + (AppSettings::preySize * 2.0f * 200.0f));

	// speed and heading
	heading = glm::vec2(.0f, -1.0f);
	speed = randomFloat(AppSettings::minPredatorVelocity, AppSettings::maxPredatorVelocity);

	target = -1;
	central = -1;
	lockOn = false;
	handling = false;
	handlingTimer = 0;
}

void Predator::calculate(std::vector<Prey>& preyAnimats) {
  
	// reset accelertion to 0 each cycle
	acceleration = glm::vec2(.0f, .0f);
  
	// only update stuff if not handling
	if (!handling) {
		// neighbour data
		glm::vec2 huntVector = glm::vec2(.0f, .0f);
		currentAttackTime++;

#if (EVOL_PARAMETERS == 1) 
			if (currentAttackTime > attackPeriod) {
				target = -1;
				handling = true;
				currentAttackTime = 0;
			}

			if ((target != -1) && (!preyAnimats[target].isDead)) {
				Prey& targetPrey = preyAnimats[target];
				float distFromTarget = glm::distance(targetPrey.position, position) - AppSettings::preySize - AppSettings::predatorSize;

				if (distFromTarget < distanceForAcceleration) isNearCatch = true;
				else isNearCatch = false;
			}
#endif
    
    // if has target
		if (target != -1) {
			// check if caught anything or target out of sight
			Prey& prey = preyAnimats[target];
			float distFromPrey = glm::distance(prey.position, position) - AppSettings::preySize - AppSettings::predatorSize;
			if (distFromPrey < .0f) distFromPrey = .0f;

			// catch attempt
			if (distFromPrey < AppSettings::catchDistance && !prey.isDead) {
				
				if (CONFUSABILITY == 1) {
					// confusability
					float random = randomFloat(.0f, 1.0f);
					int noOfConfusors = 0;
					for (Prey p2 : preyAnimats) {
						float dist = glm::distance(p2.position, position) - AppSettings::preySize - AppSettings::predatorSize;
						if (dist < .0f) dist = .0f;
						if (dist < AppSettings::confusabilitySize) noOfConfusors++;
					}

					float confusability = 0;
					if (noOfConfusors > 0) confusability = 1 / (float)noOfConfusors;

					// not confused
					if (random < confusability) {
						prey.isDead = true;
						huntCount += 1;
						currentAttackTime = 0;
					}
				}
				else {
					prey.isDead = true;
					huntCount += 1;
					currentAttackTime = 0;
				}

				target = -1;
				handling = true;
			}
			else if (CONFUSABILITY == 2) {
				huntVector = glm::normalize(prey.position - position);
				std::vector<int> targetsInAttackZone;
				std::vector<float> targetsInAttackZoneDistance;

				int index = 0;
				float min_distance = 100000.0f;
				for (Prey pr : preyAnimats) {
					if (!pr.isDead) {
						glm::vec2 p_huntVector = glm::normalize(pr.position - position);
						// angle between current target direction and animat p
						float p_alpha = std::acos((p_huntVector.x * huntVector.x + p_huntVector.y * huntVector.y));
						if (p_alpha < 0) p_alpha += 2 * std::_Pi;
						if (p_alpha < attackZoneAngle / 2) {
							float p_dist = glm::distance(pr.position, position) - AppSettings::preySize - AppSettings::predatorSize;
							if (p_dist < min_distance) min_distance = p_dist;
							targetsInAttackZone.push_back(index);
							targetsInAttackZoneDistance.push_back(p_dist);
						}
					}
					index++;
				}

				int n = targetsInAttackZone.size();
				std::vector<int> targetsInAttackZone2;
				if (n > 0) {
					for (int i = 0; i < n; i++){
						// izberemo samo animate, ki so znotraj tega polja
						if ((min_distance + 50.0f) > targetsInAttackZoneDistance[i]) {
							targetsInAttackZone2.push_back(targetsInAttackZone[i]);
						}
					}
				}
				else {
					target = -1;					currentAttackTime = 0;
				}

				n = targetsInAttackZone2.size();
				if (n > 0) target = targetsInAttackZone2[std::rand() % (n)];
				
			}
			// target out of range
			else if (distFromPrey > AppSettings::huntSize){
				target = -1;
				currentAttackTime = 0;
			}			
    }
    // else find target
    else
    {
      // find centre of nearest flock
      if (central == -1)
      {
        // copy and sort the prey by distance from predator
        std::vector<const Prey*> preyView(preyAnimats.size());
        std::transform(preyAnimats.begin(), preyAnimats.end(), preyView.begin(), [](Prey const& p) { return &p; });
        std::sort(preyView.begin(), preyView.end(), [](Prey const* a, Prey const* b)
          { return a->distanceFromPredator2 < b->distanceFromPredator2; }
        );    

        // get nearest and then the most central in his vicinity
        glm::vec2 positionOfNearest;
        bool first = false;
        float minPeripherality;
        for(auto pp = preyView.begin(); pp != preyView.end(); ++pp)
        {
          Prey const& prey = **pp;
					float distFromPred = glm::distance(prey.position, position) - AppSettings::preySize - AppSettings::predatorSize;
					if (distFromPred < .0f)
						distFromPred = .0f;
          // find nearest
          if (!first)
          {
            if (distFromPred < AppSettings::huntSize && !prey.isDead && isOnFrontSideOfSchool(prey))
            {
              first = true;
              positionOfNearest = prey.position;
              central = prey.id;
              minPeripherality = prey.peripherality;
            }
          }
          // then find most central in its vicinity
          else
          {
						float distFromNearest = glm::distance(positionOfNearest, prey.position) - AppSettings::preySize - AppSettings::predatorSize;
						if (distFromNearest < .0f)
							distFromNearest = .0f;
            // if closer than hunting zone and in nearest vicinity and smaller peripherality
            if (distFromNearest < AppSettings::cohesionSize && distFromPred < AppSettings::huntSize && 
              !prey.isDead && prey.peripherality < minPeripherality && isOnFrontSideOfSchool(prey))
            {
              minPeripherality = prey.peripherality;
              central = prey.id;
            }
          }
        }
      }

      // arppoach centre of nearest flock
      if (central != -1)
      {
        Prey& prey = preyAnimats[central];
				float distFromPrey = glm::distance(prey.position, position) - AppSettings::preySize - AppSettings::predatorSize;
				if (distFromPrey < .0f)
					distFromPrey = .0f;
        // if close enough lockOn else approach
        if (distFromPrey < lockOnDistance)
          lockOn = true;
        else
        {
          huntVector = glm::normalize(prey.position - position);
          acceleration = huntVector * AppSettings::maxPredatorForce;
        }
      }

      // lock on
      if (lockOn)
      {
        // select most peripheral in lockOnRadius
        bool targetFound = false;
        float maxPer = -std::numeric_limits<float>::max();
        for(auto p = preyAnimats.begin(); p != preyAnimats.end(); ++p)
        {
					float distFromPred = glm::distance(p->position, position) - AppSettings::preySize - AppSettings::predatorSize;
					if (distFromPred < .0f)
						distFromPred = .0f;
          if (!p->isDead && p->peripherality > maxPer && distFromPred < lockOnRadius)
          {
            maxPer = p->peripherality;
            target = p->id;
            targetFound = true;
          }
        }
        // if no possible target in lockOnRadius, cooldown (handling) and create a new attempt
        if (!targetFound)
          handling = true;
      }
    }

    // normalize and multiply with weight
    if (target != -1)
    {
      Prey pr = preyAnimats[target];
      huntVector = glm::normalize(pr.position - position);
      acceleration = huntVector * AppSettings::maxPredatorForce;
    }
  }

  // handling
  if (handling)
  {
    handlingTimer++;
    if (handlingTimer > AppSettings::handlingTime)
    {
      handling = false;
      lockOn = false;
      central = -1;
      handlingTimer = 0;
    }
  }
}

void Predator::update()
{
  // update speed and heading
  glm::vec2 velocity = getVelocity();

  if (ENERGY > 0) {
	  // added acceleration if fish still has energy
	  if (!((energy > 0) && (!isExhausted))) {
		  isExhausted = true;
		  speed = AppSettings::minPredatorVelocity; // minimum speed
		  velocity = getVelocity();
	  }

	  // if not near catch swimming x times faster than slowest speed
	  else if (!isNearCatch && EVOL_PARAMETERS == 1) {
		  speed = AppSettings::minPredatorVelocity * velocityMultiplier;
		  velocity = getVelocity();
	  }

	  else {
		  velocity += acceleration;
	  }
  }
  else {
	  velocity += acceleration;
  }


  speed = 0.0f;
  float speed2 = glm::length2(velocity);
  if (speed2 > 0.0f)
  {
    speed = std::sqrt(speed2);
    heading = velocity / speed;
  }

  // limit speed
  speed = glm::clamp(speed, AppSettings::minPredatorVelocity, AppSettings::maxPredatorVelocity);

  // if we're using oxygen consimption
  if (ENERGY == 1) {
	  oxygenConsumption = pow(62.9, (0.21 * (speed / AppSettings::predatorSize))) / 36; //mgO2 / kg h
	  oxygenConsumption = oxygenConsumption / 360;

	  if (speed > AppSettings::minPredatorVelocity) {
		  energy -= oxygenConsumption;
	  }
	  else {
		  energy += regenerationGain;
	  }
  }
  // if we're using Zheng's method
  else if (ENERGY == 2) {
	  float pi = 3.14159265358979f;
	  prevAngle_t = angle_t;
	  angle_t = atan(heading.y / heading.x);
	  float average_speed = (AppSettings::minPredatorVelocity + AppSettings::maxPredatorVelocity) / 2;
	  float energyChange = 0.001f * (1 / average_speed) * (AppSettings::minPredatorVelocity - speed) - 0.002f * (1 / pow(pi, 2)) * pow((angle_t - prevAngle_t), 2);

	  if (speed > AppSettings::minPredatorVelocity) {
		  energy += energyChange;
	  }
	  else {
		  energy += regenerationGain;
	  }
  }

  // half regenerated
  if (isExhausted && (energy >= 0.5)) {
	  isExhausted = false;
  }

  position += getVelocity();
}

bool Predator::isOnFrontSideOfSchool(Prey const& prey)
{
  if (prey.peripherality == std::numeric_limits<float>::max()) return true;
  glm::vec2 a = glm::normalize(prey.position - position);
  return glm::dot(prey.peripheralityDir, a) > AppSettings::oppositeThreshold;
}