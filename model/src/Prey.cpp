#include <cmath>
#include <limits>
#include "AppSettings.h"
#include "Predator.h"
#include "Prey.h"

Prey::Prey(int animatID)
{
  // position
  float x = randomFloat(AppSettings::screenWidth/2 - AppSettings::worldSize, 
    AppSettings::screenWidth/2 + AppSettings::worldSize);
  float y = randomFloat(AppSettings::screenHeight/2 - AppSettings::worldSize, 
    AppSettings::screenHeight/2 + AppSettings::worldSize);
  position = glm::vec2(x, y);

  // speed and heading
  x = randomFloat(-.1f * AppSettings::minPreyVelocity, .1f * AppSettings::minPreyVelocity);
  heading = glm::normalize(glm::vec2(x, -1.0f));
  speed = randomFloat(AppSettings::minPreyVelocity, AppSettings::maxPreyVelocity);

  id = animatID;
  isDead = false;
  distanceFromPredator2 = std::numeric_limits<float>::max();

  escapeDistance = 50;

  acceleration = glm::vec2(.0f, .0f);
  peripheralityDir = glm::vec2(.0f, .0f);
  peripherality = std::numeric_limits<float>::max();

  energy = 1;
  angle_t = 0;
  selfishEscape = false;
}

void Prey::calculate(Predator const& predator)
{
  // reset accelertion to 0 each cycle
  acceleration = glm::vec2(.0f, .0f);

  // neighbour data
  int escapeCount = 0;
  glm::vec2 escapeVector = glm::vec2(.0f, .0f);
  int separationCount = 0;
  glm::vec2 separationVector = glm::vec2(.0f, .0f);
  int alignmentCount = 0;
  glm::vec2 alignmentVector = glm::vec2(.0f, .0f);
  int cohesionCount = 0;
  glm::vec2 cohesionVector = glm::vec2(.0f, .0f);

  // reset peripherality
  peripheralityDir = glm::vec2(.0f, .0f);
  peripherality = std::numeric_limits<float>::max();
  int peripheralityCount = 0;

  // escape drive
  glm::vec2 ofs = position - predator.position;
  float dist2 = glm::length2(ofs);
  distanceFromPredator2 = dist2;
  if ((dist2 > 0.0) && (dist2 < AppSettings::escapeSize * AppSettings::escapeSize))
  {
		float dist = std::sqrt(dist2) - AppSettings::preySize - AppSettings::predatorSize;
		if (dist < .0f)
			dist = .0f;
		if (dist < AppSettings::escapeSize)
		{
			glm::vec2 dir = ofs / dist;
			if (!isInBlindSpot(dir))
			{
				escapeCount++;
				escapeVector = ((AppSettings::escapeSize - dist) / AppSettings::escapeSize) * dir;
			}
		}
  }

  // other drives
  for(auto ni = neighbours.begin(); ni != neighbours.end(); ++ni)
  {
    if (!isInBlindSpot(ni->dir))
    {
      if (ni->dist < AppSettings::separationSize)
      {
        separationCount++;
        glm::vec2 separation = -ni->dir;
        separation *= (AppSettings::separationSize - ni->dist) / AppSettings::separationSize;
        separationVector += separation;
      }
      else if (ni->dist < AppSettings::alignmentSize)
      {
        alignmentCount++;
        alignmentVector += ni->vel;
      }
      else if (ni->dist < AppSettings::cohesionSize)
      {
        cohesionCount++;
        cohesionVector += ni->ofs;
      }
    }

    if (ni->dist < AppSettings::peripheralitySize)
    {
      peripheralityDir += ni->dir;
      peripheralityCount ++;
    }
  }

  if (peripheralityCount != 0)
  {
    peripheralityDir /= (float)peripheralityCount;
    peripherality = glm::length(peripheralityDir);
    peripheralityDir /= peripherality;
  }

  // calculate drives
  if (escapeCount > 0)
    acceleration += (escapeVector / (float)escapeCount) * AppSettings::escapeWeight;
  if (separationCount > 0 && !selfishEscape)
    acceleration += (separationVector / (float)separationCount) * AppSettings::separationWeight;
  if (alignmentCount > 0 && !selfishEscape)
    acceleration += ((alignmentVector / (float)alignmentCount) - getVelocity()) * AppSettings::alignmentWeight;
  if (cohesionCount > 0)
    acceleration += ((cohesionVector / (float)cohesionCount)) * AppSettings::cohesionWeight;

  // clamp
  acceleration = limit(acceleration, AppSettings::maxPreyForce);

  neighbours.clear();
}

void Prey::update(Predator const& predator, std::vector<Prey> &preyAnimats)
{
  // update speed and heading
  glm::vec2 velocity = getVelocity();
  
  if (SELFISH_ESCAPE == 1 && predator.target == id) {
	  float dist = sqrt(pow(position.x - predator.position.x, 2) + pow(position.y - predator.position.y, 2));
	  if (dist < escapeDistance) selfishEscape = true;
	  else selfishEscape = false;
  }
  else selfishEscape = false;

  if (PREY_ENERGY > 0) {
	  if ((energy > 0) && (!isExhausted)) {
		  velocity += acceleration;
	  }
	  else {
		  isExhausted = true;
		  speed = AppSettings::minPreyVelocity;
		  velocity = getVelocity();
	  }
  }
  else {
	  velocity += acceleration;
  }

  speed = 0.0f;
  float speed2 = glm::length2(velocity);
  if (!selfishEscape) {
	  if (speed2 > .0f) {
		  speed = std::sqrt(speed2);
		  heading = velocity / speed;
	  }
	  // limit speed
	  speed = glm::clamp(speed, AppSettings::minPreyVelocity, AppSettings::maxPreyVelocity);
  }
  else {
	  // selfish escape
	  speed = AppSettings::maxPreyVelocity;
	  heading = glm::vec2(predator.heading.y, -predator.heading.x);
  }

  // limit speed
  speed = glm::clamp(speed, AppSettings::minPreyVelocity, AppSettings::maxPreyVelocity);

  if (PREY_ENERGY == 1) {

	  float oxygenConsumption = pow(62.9, (0.21 * ((speed / 1.5) / AppSettings::preySize))) / 36; //mgO2 / kg h
	  oxygenConsumption = oxygenConsumption / 360;

	  if (speed - glm::length2(Ui) > AppSettings::minPreyVelocity) {
		  energy -= oxygenConsumption;
		  if (energy > 1) energy = 1.0f;
	  }
	  else {
		  energy += regenerationGain;
		  if (energy > 1) energy = 1.0f;
	  }
  }

  // if we're using Zheng's method
  else if (PREY_ENERGY == 2) {

	  prevAngle_t = angle_t;
	  angle_t = atan(heading.y / heading.x);

	  ncoll = 0;
	  for (int i = 0; i < preyAnimats.size(); i++) {
		  if (id != preyAnimats[i].id && !preyAnimats[i].isDead){
			  float distance = sqrt(pow(position.x - preyAnimats[i].position.x, 2) + pow(position.y - preyAnimats[i].position.y, 2));
			  if (distance < AppSettings::preySize) ncoll += 1;
		  }
	  }

	  float average_speed = (AppSettings::minPreyVelocity + AppSettings::maxPreyVelocity) / 2;
	  float energyChange = 0.001f * (1 / average_speed) * (AppSettings::minPreyVelocity - speed)
		  - 0.002f * (1 / pow(std::_Pi, 2)) * pow((angle_t - prevAngle_t), 2)
		  - 0.002f * 1 * ncoll * std::exp(ncoll / 2);
	  if (std::abs(energyChange) > 0.9) energyChange = 0.0f;

	  if (speed - glm::length2(Ui) > AppSettings::minPreyVelocity) {
		  energy += energyChange;
	  }
	  else {
		  energy += regenerationGain;
		  if (energy > 1) energy = 1.0f;
	  }
  }

  // half regenerated
  if (isExhausted && (energy >= 0.5)) {
	  isExhausted = false;
  }

  position += getVelocity();

  if (HYDRO == 1) {
	  glm::vec2 Ui = glm::vec2(.0f, .0f);
	  for (int i = 0; i < preyAnimats.size(); i++) {
		  if (id != preyAnimats[i].id && !preyAnimats[i].isDead){
			  glm::vec2 e_jp = glm::normalize(glm::vec2(position.x - preyAnimats[i].position.x, position.y - preyAnimats[i].position.y));
			  glm::vec2 e_j0 = glm::normalize(glm::vec2(e_jp.y, -e_jp.x));
			  float p = sqrt(pow(position.x - preyAnimats[i].position.x, 2) + pow(position.y - preyAnimats[i].position.y, 2));
			  float theta = atan2(heading.y, heading.x) - atan2(preyAnimats[i].heading.y, preyAnimats[i].heading.x);

			  // polar
			  glm::vec2 e_jp_polar = glm::vec2(1, atan(e_jp.y / e_jp.x));
			  glm::vec2 e_j0_polar = glm::vec2(1, atan(e_j0.y / e_j0.x));

			  //glm::vec2 u_ji = (float)pow(10, -2) * (e_j0 * sin(theta) + e_jp * cos(theta)) / (float)(std::_Pi * pow(p, 2));
			  glm::vec2 u_ji_polar = (float)pow(10, -2) * (e_j0_polar * sin(theta) + e_jp_polar * cos(theta)) / (float)(std::_Pi * pow(p, 2));
			  //Ui += u_ji;
			  Ui += u_ji_polar;
		  }
	  }
	  position += Ui;
  }
}

bool Prey::isInBlindSpot(glm::vec2 const& animatDirection)
{
  return glm::dot(heading, animatDirection) < AppSettings::cosPreyBlindThreshold;
}