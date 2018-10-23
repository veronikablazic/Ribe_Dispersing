#pragma once
#include <vector>
#include "glm.hpp"

class Prey;
class Predator
{
  public:
    Predator();
    Predator(int animatID);
	Predator(int animatID, float parentLockOnDistance, float parentlockOnRadius, float dA, float vM, int aP);
    void calculate(std::vector<Prey>& preyAnimats);
    void update();
		void reset();
    bool isOnFrontSideOfSchool(Prey const& prey);
    glm::vec2 getVelocity() const { return speed * heading; }

    int id;
    glm::vec2 position;
    glm::vec2 heading;
    float speed;
    glm::vec2 acceleration;

    int target;
    int central;

    bool lockOn;
    float lockOnDistance;
    float lockOnRadius;

    int huntCount;
    bool handling;
    int handlingTimer;

	// new

	float regenerationGain;
	float energy;
	float angle_t;						//moving direction arctan(vy/vx) 
	float prevAngle_t;					//previous moving direction arctan(vy/vx) 
	float Ncoll;						//number of collisions
	float rFish;						//distance to other fish - for collision detection 0.15BL
	float averageSpeed;
	int averageCounter;

	float iterationsToGetRecovered;		//how many frames to get fully recovered
	float sustainedSpeed;				//sustained speed, at this speed fish can swim infinitely and can regenerate 
	float prolongedSpeed;				//a bit faster than sustained, fish will eventualy get fatigue
	float burstSpeed;					//burst speed is maximum swimming speed, fish can't swim long time at this speed
	float oxygenConsumption;			//average oxygen consumption
	float predatorDistnce;				//predator distance when fish stars to accelerating
	bool  isPeripheral;					//peripheral fish consume 10% more oxygen
	bool  isExhausted;					//if fish was exhausted it has to rest 
	float perifcount;

	float catchDistance;
	float distFromTarget;

	//speed evolution parameters
	float distanceForAcceleration;		// distance from prey to which predator is using its minimum speed
	float velocityMultiplier;			// full speed means more energy expense, why use full speed if not necessary, used instead of force multiplier
	double attackPeriod;				// max attack period in seconds
	double currentAttackTime;
	bool isNearCatch = false;
	float attackZoneAngle;

	//prey selection parameters
	double angleAttack;
	double depthAttackZone;
	double timePeriod;			//in seconds
	double prevPeriod;

	float minPreyRisk = std::numeric_limits<float>::max();
	float maxPreyRisk = std::numeric_limits<float>::lowest();
};