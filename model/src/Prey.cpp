#include <cmath>
#include <limits>
#include "AppSettings.h"
#include "Predator.h"
#include "Prey.h"
#include <algorithm>

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

  selfish = true;
  selfishTimeCount = selfishTime;

  // random
  wb = randomFloat(0.0f, 1.0f);
  wr = randomFloat(0.0f, 1.0f);
  we = randomFloat(0.0f, 1.0f);
  cn0 = randomFloat(AppSettings::minPreyVelocity, AppSettings::maxPreyVelocity);
  cn2 = randomFloat(AppSettings::minPreyVelocity, AppSettings::maxPreyVelocity);
  cn1 = randomFloat(0.0f, 50.0f);
  cn3 = randomFloat(0.0f, 50.0f);
  cr0 = randomFloat(0.0f, energy);
  cr1 = randomFloat(0.0f, 50.0f);
  selfishTime = randomInt(1, 600);
  selfishTime2 = randomInt(0, 600);
  selfishEscapeDistance = randomFloat(0.0f, AppSettings::huntSize);
  selfishProbability = randomFloat(0.0f, 1.0f);
  dot_pow = randomFloat(0.0f, 10.0f);
  dotTreshold = randomFloat(0.0f, 1.0f);
  v_b = randomFloat(AppSettings::cruisingPrey, AppSettings::maxPreyVelocity);
  v_r = randomFloat(AppSettings::minPreyVelocity, AppSettings::cruisingPrey);

}

Prey::Prey(int animatID, float wb1, float wr1, float we1, float cn01, float cn11, float cn21, float cn31, float cr01, float cr11, int st, int st2, float se, float sp, float dp, float dt, float vb, float vr)
{
	// position
	float x = randomFloat(AppSettings::screenWidth / 2 - AppSettings::worldSize, AppSettings::screenWidth / 2 + AppSettings::worldSize);
	float y = randomFloat(AppSettings::screenHeight / 2 - AppSettings::worldSize, AppSettings::screenHeight / 2 + AppSettings::worldSize);
	position = glm::vec2(x, y);

	// speed and heading
	x = randomFloat(-.1f * AppSettings::minPreyVelocity, .1f * AppSettings::minPreyVelocity);
	heading = glm::normalize(glm::vec2(x, -1.0f));
	speed = randomFloat(AppSettings::minPreyVelocity, AppSettings::maxPreyVelocity);

	id = animatID;
	isDead = false;
	//isTarget = false;

	acceleration = glm::vec2(.0f, .0f);
	peripheralityDir = glm::vec2(.0f, .0f);
	peripherality = std::numeric_limits<float>::max();

	energy = 1.0f;
	angle_t = 0;
	selfishEscape = false;
	escapeDistance = 50;

	selfish = true;
	selfishTimeCount = selfishTime;

	wb = wb1;
	wr = wr1;
	we = we1;
	cn0 = cn01;
	cn2 = cn21;
	cn1 = cn11;
	cn3 = cn31;
	cr0 = cr01;
	cr1 = cr11;
	selfishTime = st;
	selfishTime2 = st2;
	selfishEscapeDistance = se;
	selfishProbability = sp;
	dot_pow = dp;
	dotTreshold = dt;
	v_b = vb;
	v_r = vr;
}

void Prey::reset() {
	// acceleration
	acceleration = glm::vec2(.0f, .0f);

	// position
	float x = randomFloat(AppSettings::screenWidth / 2 - AppSettings::worldSize, AppSettings::screenWidth / 2 + AppSettings::worldSize);
	float y = randomFloat(AppSettings::screenHeight / 2 - AppSettings::worldSize, AppSettings::screenHeight / 2 + AppSettings::worldSize);
	position = glm::vec2(x, y);

	// speed and heading
	x = randomFloat(-.1f * AppSettings::minPreyVelocity, .1f * AppSettings::minPreyVelocity);
	heading = glm::normalize(glm::vec2(x, -1.0f));
	speed = randomFloat(AppSettings::minPreyVelocity, AppSettings::maxPreyVelocity);

	//isDead = false;
	//isTarget = false;
	//energy = 1.0f;

}

void Prey::calculate(Predator const& predator)
{

#if (SELFISH_ESCAPE == 1)
	// ko se riba odloèi za sebièni pobeg se potem nekaj èasa ne more
	if (selfishTimeout > 0) {
		selfishTimeout--;
		if (selfishTimeout == 0) selfish = true;
	}

	// sebièni pobeg traja selfishTime iteracij
	if (selfishTimeCount < selfishTime) {
		selfishTimeCount++;
		if (selfishTimeCount == selfishTime) {
			selfishEscape = false;
			selfishDirection = 0;
		}
	}
#endif


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

#if(PREY_ENERGY_PARAMS == 1) 
  //  burst -> cruise
  if (speed > cn0)
	  acceleration += wb * pow(std::min(std::max((speed - cn0) / (AppSettings::maxPreyVelocity - cn0), 0.0f), 1.0f), cn1) * glm::normalize(-getVelocity());
  // rest -> crusie
  if (speed < cn2){
	  acceleration += wr * pow(std::min(std::max((cn2 - speed) / (cn2 - AppSettings::minPreyVelocity), 0.0f), 1.0f), cn3) * glm::normalize(getVelocity());
  }
  if (energy < cr0)
	  acceleration += we * std::min(std::max((cr0 - energy / cr0), 0.0f), cr1) * glm::normalize(-getVelocity());

#endif

  // clamp
  acceleration = limit(acceleration, AppSettings::maxPreyForce);

  neighbours.clear();
}

void Prey::update(Predator const& predator, std::vector<Prey> &preyAnimats)
{
  // update speed and heading
  glm::vec2 velocity = getVelocity();
  
if (SELFISH_ESCAPE == 1) {

	  float dist = glm::distance(position, predator.position);
	  float random = randomFloat(.0f, 1.0f);

	  float dot = glm::dot(glm::normalize(predator.heading), glm::normalize(position - predator.position));
	  random = randomFloat(0.0f, 1.0f);

	  if (dist < selfishEscapeDistance && pow(dot, dot_pow) > dotTreshold) {
		  // vsaka riba samo 1x selfish pobeg
		  if (selfish) {

			  selfish = false;
			  selfishTimeout = 100;
			  selfishTimeCount = 0;

			  if (!isExhausted && random < selfishProbability) selfishEscape = true;
		  }
	  }
	  else if (selfishTimeCount == 100) {
		  selfishEscape = false;
		  selfishDirection = 0;
	  }
  }
  else selfishEscape = false;

  if (PREY_ENERGY > 0) {
	  if ((energy > 0) && (!isExhausted)) {
		  velocity += acceleration;
	  }
	  else {
		  isExhausted = true;
		  selfishEscape = false;
		  if (PREY_ENERGY_PARAMS) speed = v_r;
		  else speed = AppSettings::minPreyVelocity;
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
	  if (PREY_ENERGY_PARAMS == 1) {
		  if (!isExhausted) speed = glm::clamp(speed, v_r, v_b);
		  else speed = glm::clamp(speed, AppSettings::minPreyVelocity, v_r);
	  }
	  else speed = glm::clamp(speed, AppSettings::minPreyVelocity, AppSettings::maxPreyVelocity);
  }
  else {
	  // selfish escape
	  if (energy > 0) {
		  if (PREY_ENERGY_PARAMS == 1) speed = v_b;
		  else speed = AppSettings::maxPreyVelocity;
	  }

	  if (selfishDirection == 0) {

		  // if to the left of heading turn left, else right. if cross = 0, then random.
		  glm::vec3 ba = glm::vec3(predator.heading.x, predator.heading.y, 0);
		  glm::vec3 ca = glm::vec3(position.x - predator.position.x, position.y - predator.position.y, 0);
		  glm::vec3 cross = glm::cross(ba, ca);

		  if (cross.z > 0) {
			  heading = glm::vec2(-predator.heading.y, predator.heading.x);
			  selfishDirection = 2;
		  }
		  else if (cross.z < 0) {
			  heading = glm::vec2(predator.heading.y, -predator.heading.x);
			  selfishDirection = 1;
		  }
		  else {
			  float random = randomFloat(0, 1.0f);
			  if (random > 0.5f){
				  heading = glm::vec2(predator.heading.y, -predator.heading.x);
				  selfishDirection = 1;
			  }
			  else {
				  heading = glm::vec2(-predator.heading.y, predator.heading.x);
				  selfishDirection = 2;
			  }
		  }
	  }
	  else if (selfishDirection == 1) heading = glm::vec2(predator.heading.y, -predator.heading.x);
	  else heading = glm::vec2(-predator.heading.y, predator.heading.x);
  }

  // limit speed
  if (PREY_ENERGY_PARAMS) {
	  if (!isExhausted) speed = glm::clamp(speed, v_r, v_b);
	  else speed = glm::clamp(speed, AppSettings::minPreyVelocity, v_r);
  }
  else speed = glm::clamp(speed, AppSettings::minPreyVelocity, AppSettings::maxPreyVelocity);

#if (PREY_ENERGY == 1) 

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
#endif

  // if we're using Zheng's method
#if (PREY_ENERGY == 2) 

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

	  if (speed > AppSettings::cruisingPrey) {
		  energy += energyChange;
	  }
	  else {
		  float speedPercentage = 1 - ((speed - AppSettings::minPreyVelocity) / (AppSettings::cruisingPrey - AppSettings::minPreyVelocity));
		  energy += regenerationGain * speedPercentage;
		  if (energy > 1) energy = 1.0f;
	  }
#endif

  // half regenerated
  if (isExhausted && (energy >= 0.5)) {
	  isExhausted = false;
  }

  position += getVelocity();

#if (HYDRO == 1) 
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
#endif
}

bool Prey::isInBlindSpot(glm::vec2 const& animatDirection)
{
  return glm::dot(heading, animatDirection) < AppSettings::cosPreyBlindThreshold;
}