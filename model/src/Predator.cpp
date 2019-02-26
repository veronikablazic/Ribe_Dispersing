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
  wanderingTime = AppSettings::handlingTime;

  // new 
  energy = 1;
  regenerationGain = 30.0f / 100000.0f; //regeneration in 30'

  attackZoneAngle = (std::_Pi / 6);
  angle_t = 0;

  distFromTarget = 1000000;

  distanceForAcceleration = 765;
  attackPeriod = 248;

  wandering = false;

  currentAttackTime = 0;
  step = 0;

  if (PREY_EVOLUTION == 0) {

	  // random lock-on distances
	  lockOnDistance = randomFloat(1.0f, AppSettings::huntSize);
	  lockOnRadius = randomFloat(1.0f, AppSettings::huntSize);

	  // za part 3, ko ima plenilec energijo
	  if (PRED_ENERGY_PARAMS == 1) {
		  wb = randomFloat(0.0f, 1.0f);
		  wr = randomFloat(0.0f, 1.0f);
		  we = randomFloat(0.0f, 1.0f);
		  cn0 = randomFloat(AppSettings::minPredatorVelocity, AppSettings::maxPredatorVelocity);
		  cn2 = randomFloat(AppSettings::minPredatorVelocity, AppSettings::maxPredatorVelocity);
		  cn1 = randomFloat(0.0f, 50.0f);
		  cn3 = randomFloat(0.0f, 50.0f);
		  cr0 = randomFloat(0.0f, energy);
		  cr1 = randomFloat(0.0f, 50.0f);

		  attackPeriod = randomInt(1, AppSettings::noOfSteps);
		  wanderingTime = randomInt(0, AppSettings::noOfSteps);
		  //velocityMultiplier = randomFloat(1.0f, 3.0f);
		  distanceForAcceleration = randomFloat(1.0f, AppSettings::huntSize);
		  //nextDecision = restPeriod;

		  v_b = randomFloat(AppSettings::cruisingSpeed, AppSettings::maxPredatorVelocity);
		  v_r = randomFloat(AppSettings::minPredatorVelocity, AppSettings::cruisingSpeed);
	  }
  }
  // v part 2 ima vse parametre fiksno dolocene, te 3 k jih ma
  else {

	  // for conf = 0 and conf = 2

	  lockOnDistance = randomFloat(1.0f, AppSettings::huntSize);
	  lockOnRadius = randomFloat(1.0f, AppSettings::huntSize);

	  // end for conf = 0
  }
}

Predator::Predator(int animatID, float parentLockOnDistance, float parentlockOnRadius, float dA, int aP, int rP, float wb1, float wr1, float we1, float cn01, float cn11, float cn21, float cn31, float cr01, float cr11, float vb, float vr)
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
  attackPeriod = aP;
  currentAttackTime = 0;

  wb = wb1;
  wr = wr1;
  we = we1;
  cn0 = cn01;
  cn2 = cn21;
  cn1 = cn11;
  cn3 = cn31;
  cr0 = cr01;
  cr1 = cr11;

  wanderingTime = rP;

  v_b = vb;
  v_r = vr;
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
	wanderingTimer = 0;

	energy = 1.0f;
}

void Predator::calculate(std::vector<Prey>& preyAnimats) {
  
	// reset accelertion to 0 each cycle
	acceleration = glm::vec2(.0f, .0f);
  
	// only update stuff if not handling
	if (!(handling || wandering)) {
		// neighbour data
		glm::vec2 huntVector = glm::vec2(.0f, .0f);
		currentAttackTime++;

		#if (PRED_ENERGY_PARAMS == 1) 
					if (currentAttackTime > attackPeriod) {
						target = -1;
						wandering = true;
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
				
				#if (CONFUSABILITY == 1) 
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
							handling = true;
						}
						else {
							wandering = true;
						}
				#endif
				#if (CONFUSABILITY != 1) 
						prey.isDead = true;
						huntCount += 1;
						currentAttackTime = 0;
						handling = true;
				#endif

				target = -1;
				isNearCatch = false;
			}

			#if (CONFUSABILITY == 2)
				huntVector = glm::normalize(prey.position - position);
				std::vector<int> targetsInAttackZone;

				int index = 0;
				glm::vec2 p_huntVector;
				float p_cos_alpha;
				float p_dist;

				for (Prey pr : preyAnimats) {
					if (!pr.isDead) {
						p_huntVector = glm::normalize(pr.position - position);
						p_cos_alpha = glm::dot(p_huntVector, huntVector);

						if (p_cos_alpha > cos(attackZoneAngle / 2.0f)) {
							p_dist = glm::distance(prey.position, pr.position) - AppSettings::preySize * 2;
							if (p_dist < AppSettings::confusabilitySize / 2) {
								targetsInAttackZone.push_back(index);
							}
						}
					}
					index++;
				}

				int n = targetsInAttackZone.size();
				if (n > 0) {
					target = targetsInAttackZone[std::rand() % (n)];
				}			
		#endif
			// target out of range
		if (distFromPrey > AppSettings::huntSize){
			target = -1;
			currentAttackTime = 0;
			if (PRED_ENERGY_PARAMS == 1) wandering = true;
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
    if (target != -1) {
      Prey pr = preyAnimats[target];
      huntVector = glm::normalize(pr.position - position);
	  if (!isExhausted) {

		  // for part 1 we want to ignore energy consumption, therefore maxForce
		  if (PRED_ENERGY_PARAMS == 0) {
			  acceleration = huntVector * AppSettings::maxPredatorForce;
		  }
		  // else we take the base huntVector and add the new drives
		  else {
			  acceleration = huntVector;
			  //  burst -> cruise
			  if (speed > cn0)
				  acceleration += pow(std::min(std::max((speed - cn0) / (AppSettings::maxPredatorVelocity - cn0), 0.0f), 1.0f), cn1) * glm::normalize(-getVelocity());
			  // rest -> crusie
			  if (speed < cn2){
				  float wn = pow(std::min(std::max((cn2 - speed) / (cn2 - AppSettings::minPredatorVelocity), 0.0f), 1.0f), cn3);
				  acceleration += wn * glm::normalize(getVelocity());
			  }
			  if (energy < cr0){
				  acceleration += we * std::min(std::max((cr0 - energy / cr0), 0.0f), cr1) * glm::normalize(-getVelocity());
			  }
		  }
	  }
    }
  }

  // handling
  if (handling || wandering)
  {

	  isNearCatch = false;

	  handlingTimer++;
	  if (handlingTimer > AppSettings::handlingTime) {
		  handling = false;
		  handlingTimer = 0;
		  target = -1;
		  currentAttackTime = 0;
	  }

	  wanderingTimer++;
	  if (wanderingTimer > wanderingTime) {
		  wandering = false;
		  wanderingTimer = 0;
		  target = -1;
		  currentAttackTime = 0;
	  }

	  nextDecision--;
	  if (nextDecision == 0) {

		  glm::vec2 direction = glm::vec2(randomFloat(-100.0f, 100.0f), randomFloat(-100.0f, 100.0f));
		  glm::vec2 unit_temp = glm::normalize(direction);
		  nextDecision = randomInt(10, 30);

		  float p_cos_alpha = glm::dot(heading, unit_temp);

		  // je znotraj pahljace
		  if (p_cos_alpha > cos(std::_Pi / 2)) {
			  unit = unit_temp;
		  }
		  else unit = -unit_temp;
	  }

	  if (PRED_ENERGY_PARAMS == 0) acceleration = unit * AppSettings::maxPredatorForce * 0.2f;
	  else {
			acceleration = unit;
			//  burst -> cruise
			if (speed > cn0) acceleration += pow(std::min(std::max((speed - cn0) / (AppSettings::maxPredatorVelocity - cn0), 0.0f), 1.0f), cn1) * glm::normalize(-getVelocity());
			// rest -> crusie
			if (speed < cn2){
				float wn = pow(std::min(std::max((cn2 - speed) / (cn2 - AppSettings::minPredatorVelocity), 0.0f), 1.0f), cn3);
				acceleration += wn * glm::normalize(getVelocity());
			}
			if (energy < cr0) acceleration += we * std::min(std::max((cr0 - energy / cr0), 0.0f), cr1) * glm::normalize(-getVelocity());
	  }
  }
}

void Predator::update(std::vector<Prey>& preyAnimats)
{
  // update speed and heading
  glm::vec2 velocity = getVelocity();

  if (ENERGY > 0) {
	  if (isNearCatch && (PRED_ENERGY_PARAMS == 1)) speed = v_b;
	  else if (wandering) speed = v_r;
	  else speed = AppSettings::cruisingSpeed;
  }

  if (PRED_ENERGY_PARAMS == 1){
	  // added acceleration if fish still has energy
	  if (energy <= 0) {
		  isExhausted = true;
		  speed = AppSettings::minPredatorVelocity; // minimum speed
	  }
	  else {
		  step++;
	  }
  }

  velocity = getVelocity();
  if (!isExhausted) velocity += acceleration;

  speed = 0.0f;
  float speed2 = glm::length2(velocity);
  if (speed2 > 0.0f) {
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

  if (energy > 1.0f) energy = 1.0f;
  position += getVelocity();
}

bool Predator::isOnFrontSideOfSchool(Prey const& prey)
{
  if (prey.peripherality == std::numeric_limits<float>::max()) return true;
  glm::vec2 a = glm::normalize(prey.position - position);
  return glm::dot(prey.peripheralityDir, a) > AppSettings::oppositeThreshold;
}