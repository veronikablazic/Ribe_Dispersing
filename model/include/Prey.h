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
	Prey(int animatID, float wb1, float wr1, float we1, float cn01, float cn11, float cn21, float cn31, float cr01, float cr11, int st, int st2, float se, float sp, float dp, float dt, float vb, float vr);
    void calculate(Predator const& predator);
	void update(Predator const& predator, std::vector<Prey> &preyAnimats);
    bool isInBlindSpot(glm::vec2 const& animatDirection);
    glm::vec2 getVelocity() const { return speed * heading; }
	void reset();

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
	bool selfish;
	int selfishTimeout;
	int selfishTimeCount = 100;
	int selfishDirection = 0;
	glm::vec2 Ui = glm::vec2(.0f, .0f);

	float wb, wr, we;
	float cn0, cn1, cn2, cn3;
	float cr0, cr1;

	int selfishTime = 30;					// koliko iteracij je plen v selfish nacinu
	int selfishTime2 = 100;					// koliko èasa plen ne more veè preiti v selfish
	float selfishEscapeDistance = 50.0f;	// radij od plenilca, pri katerem se prozi
	float selfishProbability = 0.5f;		// 5% verjetnost v vsaki iteraci, ki se manjsa z dot.
	float dot_pow = 5.0f;					// potenca za vpliv dot produkta na P(selfish)
	float dotTreshold = 0.9f;				// meja za dot, po kateri se proži selfish

	float v_b = AppSettings::maxPreyVelocity;
	float v_r = AppSettings::minPreyVelocity;
};
