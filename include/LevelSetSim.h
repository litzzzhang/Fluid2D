#pragma once

#include <Eigen/Sparse>
#include <vector>
#include <glm/glm.hpp>

class LevelSetSim2D {
public:
	float width_;
	int ni_, nj_;
	float h_;

	// simulation data
	std::vector<float> u_, v_, u_temp_, v_temp_;
	
	// pressure solver data using eigen data container
	Eigen::VectorXd p_; // p is not real pressure, it contains const dt/rho
	Eigen::VectorXd rhs_;
	Eigen::SparseMatrix<double> A_; // to-do figure out to use column major or row major or both

	// geometry data
	std::vector<float> boundary_sdf_, liquid_sdf_;
	// weights represents how much portion is not in the solid boundary
	std::vector<float> u_weights_, v_weights_;
	// valid marks velocity which stands for liquid velocity, exclude pure air and solid boundary
	std::vector<char> u_valid_, v_valid_;

	LevelSetSim2D() = delete;
	LevelSetSim2D(float width, int nx, int ny);

	// simulation
	glm::vec2 GetVelocity(glm::vec2 pos) const;
	float GetValue(glm::vec2 pos) const;
	void Advance(float dt);
	// geometry
	void SetBoundary(float sdf_function(const glm::vec2&));
	void SetLiquid(float sdf_function(const glm::vec2&));



private:
	void AddForce(float dt);
	void Advect(float dt);
	void Project(float dt);
	void ConstrainBoundary();

	float CFL() const;
	glm::vec2 TraceRK2(const glm::vec2& pos, float dt) const;

	// helper function for pressure projection
	void ComputeSDF();
	void ComputeWeights();
	void SolvePressure();
	void ApplyPressureGradient();

	void ExtrapolateToBoundary(std::vector<float>& velocity_field, int vel_ni, int vel_nj,
		std::vector<char>& valid);
};
