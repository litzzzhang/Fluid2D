#include "LevelSetSim.h"
#include "util.h"

#include <algorithm>

LevelSetSim2D::LevelSetSim2D(float width, int ni, int nj)
: width_(width), ni_(ni), nj_(nj), h_(width / (float) ni)
{
	u_.resize((ni_ + 1) * nj_);
	u_temp_.resize((ni_ + 1) * nj_);
	u_valid_.resize((ni_ + 1) * nj_);
	v_.resize((nj_ + 1) * ni_);
	v_temp_.resize((nj_ + 1) * ni_);
	v_valid_.resize((nj_ + 1) * ni_);

	p_.resize(ni_ * nj_);
	
	boundary_sdf_.resize((ni_ + 1) * (nj_ + 1));
	liquid_sdf_.resize(ni_ * nj_);

	// initialize velocity field with zero velocity
	u_.assign((ni_ + 1) * nj_, 0.0f);
	v_.assign((nj_ + 1) * ni_, 0.0f);
}

glm::vec2 LevelSetSim2D::GetVelocity(glm::vec2 pos) const
{
	// subtract offset in order to convert MAC grid to grid centered
	float u_sample = InterpolateValueOnGrid(pos / h_ - glm::vec2(0.0f, 0.5f), u_, ni_ + 1, nj_);
	float v_sample = InterpolateValueOnGrid(pos / h_ - glm::vec2(0.5f, 0.0f), v_, ni_, nj_ + 1);

	return glm::vec2(u_sample, v_sample);
}

void LevelSetSim2D::Advance(float dt)
{
	float t = 0.0f;

	while (t < dt) {
		float substep = CFL();
		if (t + substep > dt) {
			substep = dt - t;
		}
		printf("Taking substep of size %f (to %0.3f%% of the frame)\n", substep, 100 * (t + substep) / dt);

		Advect(substep);
		AddForce(substep);
		Project(substep);

		ExtrapolateToBoundary(u_, ni_ + 1, nj_, u_valid_, ni_ + 1, nj_);
		ExtrapolateToBoundary(v_, ni_, nj_ + 1, v_valid_, ni_, nj_ + 1);

		ConstrainBoundary();
	}
}

void LevelSetSim2D::SetBoundary(float sdf_function(const glm::vec2&))
{
	for (int j = 0; j < nj_ + 1; j++) {
		for (int i = 0; i < ni_ + 1; i++) {
			boundary_sdf_[i + (ni_ + 1) * j] = sdf_function({ i * h_, j * h_ });
		}
	}

}

void LevelSetSim2D::SetLiquid(float sdf_function(const glm::vec2&))
{

	for (int j = 0; j < nj_; j++) {
		for (int i = 0; i < ni_; i++) {
			liquid_sdf_[i + j * (nj_)] = sdf_function({ (i + 0.5f) * h_, (j + 0.5f) * h_});
		}
	}
}

void LevelSetSim2D::AddForce(float dt)
{
	for (int j = 0; j < nj_ + 1; j++) {
		for (int i = 0; i < ni_; i++) {
			v_[i + j * ni_] += dt * -9.8f;
		}
	}
}

void LevelSetSim2D::Advect(float dt)
{
	u_temp_.assign((ni_ + 1) * nj_, 0.0f);
	v_temp_.assign(ni_ * (nj_ + 1), 0.0f);

	for (int j = 0; j < nj_; j++) {
		for (int i = 0; i < ni_ + 1; i++) {
			glm::vec2 pos(i * h_, (j + 0.5f) * h_); // MAC grid
			pos = TraceRK2(pos, -dt);
			u_temp_[i + (ni_ + 1) * j] = GetVelocity(pos)[0];
		}
	}

	for (int j = 0; j < nj_ + 1; j++) {
		for (int i = 0; i < ni_; i++) {
			glm::vec2 pos((i + 0.5f) * h_, j * h_); // MAC grid
			pos = TraceRK2(pos, -dt);
			v_temp_[i + ni_ * j] = GetVelocity(pos)[1];
		}
	}

	// possible trial u_ = std::move(u_temp_);
	// if the u_temp is really needed? can we just allocate resource in the advect step?
	u_.swap(u_temp_);
	v_.swap(v_temp_);
}

void LevelSetSim2D::Project(float dt)
{
	ComputeSDF();
	ComputeWeights();
	SolvePressure();
	ApplyPressureGradient(dt);
}

void LevelSetSim2D::ConstrainBoundary()
{
	u_temp_.swap(u_);
	v_temp_.swap(v_);

	for (int j = 0; j < nj_; j++) {
		for (int i = 0; i < ni_ + 1; i++) {
			// for velocity in the solid boundary
			if (u_weights_[i + (ni_ + 1) * j] == 0) {
				glm::vec2 pos(i * h_, (j + 0.5f) * h_); // MAC grid
				glm::vec2 vel = GetVelocity(pos);
				glm::vec2 surface_normal = InterpolateGradientOnGrid(pos / h_, boundary_sdf_, ni_ + 1, nj_ + 1);;
				
				if (glm::length(surface_normal) > 0.0f) {
					surface_normal = glm::normalize(surface_normal);
					float vel_n = glm::dot(surface_normal, vel);
					vel = vel - vel_n * surface_normal;
					u_temp_[i + j * (ni_ + 1)] = vel.x;
				}
				else {
					u_temp_[i + j * (ni_ + 1)] = 0.0f;
				}
			}
		}
	}

	for (int j = 0; j < nj_ + 1; j++) {
		for (int i = 0; i < ni_; i++) {
			// for velocity in the solid boundary
			if (v_weights_[i + ni_ * j] == 0) {
				glm::vec2 pos((i + 0.5f) * h_, j * h_); // MAC grid
				glm::vec2 vel = GetVelocity(pos);
				glm::vec2 surface_normal = InterpolateGradientOnGrid(pos / h_, boundary_sdf_, ni_ + 1, nj_ + 1);;

				if (glm::length(surface_normal) > 0.0f) {
					surface_normal = glm::normalize(surface_normal);
					float vel_n = glm::dot(surface_normal, vel);
					vel = vel - vel_n * surface_normal;
					v_temp_[i + j * ni_] = vel.x;
				}
				else {
					v_temp_[i + j * ni_] = 0.0f;
				}
			}
		}
	}

	u_temp_.swap(u_);
	v_temp_.swap(v_);
}

float LevelSetSim2D::CFL() const
{
	float max_velocity = 0.0f;
	for (const float& u : u_) {
		if (abs(u) > max_velocity)
			max_velocity = abs(u);
	}
	for (const float& v : v_) {
		if (abs(v) > max_velocity)
			max_velocity = abs(v);
	}

	max_velocity += sqrt(5.0f * h_ * 9.8f); // in case of zero divide
	return h_ / max_velocity;
}

glm::vec2 LevelSetSim2D::TraceRK2(const glm::vec2& pos, float dt) const
{
	glm::vec2 velocity = GetVelocity(pos);

	glm::vec2 mid_point = pos + 0.5f * dt * velocity;

	velocity = GetVelocity(mid_point);

	return pos + dt * velocity;
}

void LevelSetSim2D::ComputeSDF()
{
	std::vector<float> temp_liquid_sdf(ni_ * nj_, 100.0f);

	// iterate current liquid sdf, find the boundary(sign change)
	// update sdf in temp_liquid_sdf on boundary grid
	for (int j = 0; j < nj_; j++) {
		for (int i = 0; i < ni_; i++) {
			float center = liquid_sdf_[i + j * ni_];
			float left = i == 0 ? center : liquid_sdf_[i - 1 + j * ni_];
			float right = i == ni_ - 1 ? center : liquid_sdf_[i + 1 + j * ni_];
			float down = j == 0 ? center : liquid_sdf_[i + (j - 1) * ni_];
			float up = j == nj_ - 1 ? center : liquid_sdf_[i + (j + 1) * ni_];

			std::vector<float> neighbors{left, right, down, up};

			float dist = 100.0f;
			auto calculate_ratio = [&](float neighbor) {
					return abs(center) / (center - neighbor);
			};

			for (const float& neighbor : neighbors) {
				if (neighbor * center < 0.0f) { 
				// for neighbors who don't exist, let them be the same as center
				// so that their ratio will not be calculate
					if (calculate_ratio(neighbor) * h_ < dist) {
						dist = calculate_ratio(neighbor);
					}
				}
			}
			if (dist < 100.0f) // if we do update sdf
				temp_liquid_sdf[i + j * nj_] = dist;
		}
	}
	// loop through the grid to propagate the distance with eikonal equation
	// sort and solve eikonal equation with 1, 2, (3) respectively

	for (int iteration = 0; iteration < 2; iteration++) {//fast sweep
		for (int j = 0; j < nj_; j++) {
			for (int i = 0; i < ni_; i++) {
				float center = liquid_sdf_[i + j * ni_];
				// if grid not exist, initialize with infinite value(100.0f)
				// so that we will definitely select existing grid for Eikonal equation
				float left = i == 0 ? 100.0f : liquid_sdf_[i - 1 + j * ni_];
				float right = i == ni_ - 1 ? 100.0f : liquid_sdf_[i + 1 + j * ni_];
				float down = j == 0 ? 100.0f : liquid_sdf_[i + (j - 1) * ni_];
				float up = j == nj_ - 1 ? 100.0f : liquid_sdf_[i + (j + 1) * ni_];

				float x_min = std::min(left, right);
				float y_min = std::min(up, down);

				if (x_min * y_min > 0.0f) {
					int is_negative = (x_min < 0.0f);
					x_min = abs(x_min);
					y_min = abs(y_min);
					float dist = std::min(x_min, y_min) + h_;
					if (dist > std::max(x_min, y_min)) {
						dist = 0.5f * (x_min + y_min + std::sqrt(2 * h_ * h_ - std::pow(x_min- y_min, 2)));
					}
					if (dist < temp_liquid_sdf[i + j * ni_]) {
						temp_liquid_sdf[i + j * ni_] = pow(-1, is_negative) * dist;
					}
				}

			}

			for (int i = ni_ - 1; i > 0; i--) {
				float center = liquid_sdf_[i + j * ni_];
				// if grid not exist, initialize with infinite value(100.0f)
				// so that we will definitely select existing grid for Eikonal equation
				float left = i == 0 ? 100.0f : liquid_sdf_[i - 1 + j * ni_];
				float right = i == ni_ - 1 ? 100.0f : liquid_sdf_[i + 1 + j * ni_];
				float down = j == 0 ? 100.0f : liquid_sdf_[i + (j - 1) * ni_];
				float up = j == nj_ - 1 ? 100.0f : liquid_sdf_[i + (j + 1) * ni_];


				float x_min = std::min(left, right);
				float y_min = std::min(up, down);

				if (x_min * y_min > 0.0f) {
					int is_negative = (x_min < 0.0f);
					x_min = abs(x_min);
					y_min = abs(y_min);
					float dist = std::min(x_min, y_min) + h_;
					if (dist > std::max(x_min, y_min)) {
						dist = 0.5f * (x_min + y_min + std::sqrt(2 * h_ * h_ - std::pow(x_min- y_min, 2)));
					}
					if (dist < temp_liquid_sdf[i + j * ni_]) {
						temp_liquid_sdf[i + j * ni_] = pow(-1, is_negative) * dist;
					}
				}
			}
		}

		for (int j = nj_ - 1; j > 0; j--) {
			for (int i = 0; i < ni_; i++) {
				float center = liquid_sdf_[i + j * ni_];
				// if grid not exist, initialize with infinite value(100.0f)
				// so that we will definitely select existing grid for Eikonal equation
				float left = i == 0 ? 100.0f : liquid_sdf_[i - 1 + j * ni_];
				float right = i == ni_ - 1 ? 100.0f : liquid_sdf_[i + 1 + j * ni_];
				float down = j == 0 ? 100.0f : liquid_sdf_[i + (j - 1) * ni_];
				float up = j == nj_ - 1 ? 100.0f : liquid_sdf_[i + (j + 1) * ni_];

				float x_min = std::min(left, right);
				float y_min = std::min(up, down);

				if (x_min * y_min > 0.0f) {
					int is_negative = (x_min < 0.0f);
					x_min = abs(x_min);
					y_min = abs(y_min);
					float dist = std::min(x_min, y_min) + h_;
					if (dist > std::max(x_min, y_min)) {
						dist = 0.5f * (x_min + y_min + std::sqrt(2 * h_ * h_ - std::pow(x_min - y_min, 2)));
					}
					if (dist < temp_liquid_sdf[i + j * ni_]) {
						temp_liquid_sdf[i + j * ni_] = pow(-1, is_negative) * dist;
					}
				}

			}

			for (int i = ni_ - 1; i > 0; i--) {
				float center = liquid_sdf_[i + j * ni_];
				// if grid not exist, initialize with infinite value(100.0f)
				// so that we will definitely select existing grid for Eikonal equation
				float left = i == 0 ? 100.0f : liquid_sdf_[i - 1 + j * ni_];
				float right = i == ni_ - 1 ? 100.0f : liquid_sdf_[i + 1 + j * ni_];
				float down = j == 0 ? 100.0f : liquid_sdf_[i + (j - 1) * ni_];
				float up = j == nj_ - 1 ? 100.0f : liquid_sdf_[i + (j + 1) * ni_];

				float x_min = std::min(left, right);
				float y_min = std::min(up, down);

				if (x_min * y_min > 0.0f) {
					int is_negative = (x_min < 0.0f);
					x_min = abs(x_min);
					y_min = abs(y_min);
					float dist = std::min(x_min, y_min) + h_;
					if (dist > std::max(x_min, y_min)) {
						dist = 0.5f * (x_min + y_min + std::sqrt(2 * h_ * h_ - std::pow(x_min - y_min, 2)));
					}
					if (dist < temp_liquid_sdf[i + j * ni_]) {
						temp_liquid_sdf[i + j * ni_] = pow(-1, is_negative) * dist;
					}
				}
			}

		}
	}

	liquid_sdf_.swap(temp_liquid_sdf);
}

void LevelSetSim2D::ComputeWeights()
{
	float weight = 0.0f;
	for (int j = 0; j < nj_; j++) {
		for (int i = 0; i < ni_ + 1; i++) {
			// index will not exceed because boundary sdf is defined on the grid nodal 
			// u_weights is defined on the face center
			weight = 1 - FractionInside(boundary_sdf_[i + j * ni_],
												boundary_sdf_[i + (j + 1) * ni_]);
			u_weights_[i + j * ni_] = glm::clamp(weight, 0.0f, 1.0f);
		}
	}

	for (int j = 0; j < nj_ + 1; j++) {
		for (int i = 0; i < ni_; i++) {
			// index will not exceed because boundary sdf is defined on the grid nodal 
			// u_weights is defined on the face center
			weight = 1 - FractionInside(boundary_sdf_[i + j * ni_],
												boundary_sdf_[i + 1 + j * ni_]);
			v_weights_[i + j * ni_] = glm::clamp(weight, 0.0f, 1.0f);
		}
	}
}

void LevelSetSim2D::SolvePressure()
{
	unsigned int size = ni_ * nj_;
	unsigned int active_size = 5 * size - 2 * ni_ - 2 * nj_;
	
	if (rhs_.size() != size) {
		rhs_.resize(size);
		p_.resize(size);
		A_.resize(size, size);
		// exclude edge case
		A_.reserve(active_size);
	}
	
	std::vector<Eigen::Triplet<double>> TripletList;
	TripletList.reserve(active_size * sizeof(Eigen::Triplet<double>));
	rhs_.setZero();
	p_.setZero();
	A_.setZero();

	double invh = 1 / h_;
	double invh2 = invh * invh;
	for (int j = 0; j < nj_; j++){
		for (int i = 0; i < ni_; i++){
			int index = i + j * ni_;
			float center_sdf = liquid_sdf_[index];
			double center_value = 0.0;

			if (center_sdf < 0.0f){
				double term = 0.0;
				// right neighbor, excluding from edge case
				if (i < ni_ - 1){
					double right_value = 0.0;
					term = u_weights_[index + 1] * invh2;
					float right_sdf = liquid_sdf_[index + 1];
					if (right_sdf < 0.0f){ // if right grid is also liquid
						center_value += term;
						right_value -= term;
					}
					else{
						double theta = FractionInside(right_sdf, center_sdf);
						theta = std::max(theta, 0.01);
						center_value += term / theta;
					}
					if (right_value != 0.0)
						TripletList.emplace_back(index, index + 1, right_value);
					rhs_(index) -= u_weights_[index + 1] * u_[index + 1] * invh;
				}
				// left neighbor, excluding from edge cases
				if (i > 0){
					term = u_weights_[index - 1] * invh2;
					float left_sdf = liquid_sdf_[index - 1];
					double left_value = 0.0;
					if (left_sdf < 0.0f){
						center_value += term;
						left_value -= term;
					}else{
						double theta = FractionInside(left_sdf, center_sdf);
						theta = std::max(theta, 0.01);
						center_value += term / theta;
					}
					if (left_value != 0.0)
						TripletList.emplace_back(index, index - 1, left_value);
					rhs_(index) += u_weights_[index - 1] * u_[index - 1] * invh;
				}
				// up neighbor, excluding from edge cases
				if (j < nj_ - 1){
					term = v_weights_[index + nj_] *invh2;
					float up_sdf = liquid_sdf_[index + nj_];
					double up_value = 0.0;
					if (up_sdf < 0.0f){
						center_value += term;
						up_value -= term;
					}else{
						double theta = FractionInside(up_sdf, center_sdf);
						theta = std::max(theta, 0.01);
						center_value += term / theta;
					}
					if (up_value != 0.0)
						TripletList.emplace_back(index, index + nj_, up_value);
					rhs_(index) -= v_weights_[index + nj_] * v_[index + nj_] * invh;
				}
				// down neighbor, excluding from edge cases
				if (j > 0){
					term = v_weights_[index - nj_] * invh2;
					float down_sdf = liquid_sdf_[index - nj_];
					double down_value = 0.0;
					if (down_sdf < 0.0f){
						center_value += term;
						down_value -= term;
					}else{
						double theta = FractionInside(down_sdf, center_sdf);
						theta = std::max(theta, 0.01);
						center_value += term / theta;
					}
					if (down_value != 0.0)
						TripletList.emplace_back(index, index - nj_, down_value);
					rhs_(index) += v_weights_[index - nj_] * v_[index - nj_] * invh;
				}
				TripletList.emplace_back(index, index, center_value);
			}else{
				// for non-liquid grid, just set diag to be one so the pressure here equals 0
				TripletList.emplace_back(index, index, 1.0);
			}
		}
	}
	// construction down
	A_.setFromTriplets(TripletList.begin(), TripletList.end());

	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> cg_solver;
	cg_solver.compute(A_);
	if (cg_solver.info() != Eigen::Success){
		printf("Decomposition Failed!");	
	}
	p_ = cg_solver.solve(rhs_);
	if (cg_solver.info() != Eigen::Success){
		printf("Solve Failed!");	
	}
}

void LevelSetSim2D::ApplyPressureGradient(float dt)
{
}

void LevelSetSim2D::ExtrapolateToBoundary(std::vector<float>& velocity_field, int vel_ni, int vel_nj,
									std::vector<char>& valid, int marker_ni, int marker_nj)
{
	std::vector<float> temp_field = velocity_field;
	std::vector<char> prev_valid(marker_ni * marker_nj);

	for (int layer = 0; layer < 10; layer++) { // to-do: let CFL decide the iteration
		prev_valid = valid;

		// start from inner to get rid of access out of boundary
		for (int j = 1; j < vel_nj - 1; j++) {
			for (int i = 1; i < vel_ni - 1; i++) {
				float sum = 0.0f;
				int count = 0;

				if (!prev_valid[i + vel_ni * j]) {
					if (prev_valid[i + 1 + vel_ni * j]) {
						count++;
						sum += velocity_field[i + 1 + vel_ni * j];
					}
					if (prev_valid[i - 1 + vel_ni * j]) {
						count++;
						sum += velocity_field[i - 1 + vel_ni * j];
					}
					if (prev_valid[i + vel_ni * (j + 1)]) {
						count++;
						sum += velocity_field[i + vel_ni * (j + 1)];
					}
					if (prev_valid[i + vel_ni * (j - 1)]) {
						count++;
						sum += velocity_field[i + vel_ni * (j - 1)];
					}

					if (count > 0) {
						temp_field[i + j * ni_] = sum / (float)count;
						valid[i + j * ni_] = 1;
					}
				}
			}
		}
		velocity_field = temp_field;
	}
}
