#pragma once
#include <vector>
#include "glm.hpp"

struct neighbour_info
{
  neighbour_info(glm::vec2 const& Ofs, glm::vec2 const& Vel, float Dist)
    : ofs(Ofs), vel(Vel), dir(Ofs / Dist), dist(Dist) {}

  glm::vec2 ofs;
  glm::vec2 vel;
  glm::vec2 dir;
  float dist;
};

class Prey 
{
  public:
    Prey(int animatID);
    void calculate(Predator const& predator);
	void update(Predator const& predator, std::vector<Prey> &preyAnimats);
    bool isInBlindSpot(glm::vec2 const& animatDirection);
    glm::vec2 getVelocity() const { return speed * heading; }

    bool isDead;
    int id;
    glm::vec2 position;
    glm::vec2 heading;
    float speed;
    glm::vec2 acceleration;

    float distanceFromPredator2;

    glm::vec2 peripheralityDir;
    float peripherality;

    std::vector<neighbour_info> neighbours;

	//evolution parameter
	float escapeDistance;				//predator distance when fish stars accelerating

	float energy;
	bool isExhausted = false;
	float regenerationGain = 30.0f / 10000.0f; //regeneration in 30'
	float prevAngle_t;
	float angle_t;
	double ncoll = 0;

	bool selfishEscape = false;
	glm::vec2 Ui = glm::vec2(.0f, .0f);
};
