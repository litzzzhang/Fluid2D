#include "LevelSetSim.h"

const float r0 = 0.5f;
const glm::vec2 center = glm::vec2(0.5f);
const float timestep = 0.01f;

float boundary_sdf(const glm::vec2& pos) {
	float dist = glm::length(center - pos);
	return -(dist - r0);
}

float liquid_sdf(const glm::vec2& pos) {
	return 0.0f;
}

int main() {
	LevelSetSim2D sim(1.0f, 100, 100);

	// initialization boundary
	sim.SetBoundary(boundary_sdf);
	sim.SetLiquid(liquid_sdf);

	for (int i = 0; i < 1000; i++) {
		sim.Advance(timestep);
	}

	return 0;
}