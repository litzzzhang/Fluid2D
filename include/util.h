#pragma once
#include <cmath>
#include <glm/glm.hpp>
#include <vector>

void GetBaryCentric(float x, int& i, float& f, int i_low_bound, int i_high_bound) {
	// x is the coordinate after grid normalization
	// i represents for the left point in the linear interpolation
	// low bound and high bound represents for the minimal/maximal grid index of possible left point
	float x_floor = std::floor(x);
	i = (int)x_floor;

	if (i < i_low_bound) {
		i = i_low_bound;
		f = 0.0f;
	}
	else if (i > i_high_bound) {
		i = i_high_bound;
		f = 1.0f;
	}
	else {
		f = x - x_floor;
	}
}

float Lerp(float y0, float y1, float f) {
	return (1 - f) * y0 + f * y1;
}

float BiLerp(float y00, float y01, float y10, float y11, float fx, float fy) {
	return Lerp(Lerp(y00, y01, fy), Lerp(y10, y11, fy), fx);
}

float InterpolateValueOnGrid(const glm::vec2& pos_normalized, const std::vector<float>& field, const int ni, const int nj){
	int i, j;
	float fx, fy;

	GetBaryCentric(pos_normalized.x, i, fx, 0, ni - 2); // max index is field[ni - 1], left point max is ni - 2
	GetBaryCentric(pos_normalized.y, j, fy, 0, nj - 2); // max index is field[nj - 1], left point max is nj - 2

	return BiLerp(field[i + ni * j], field[i + ni * (j + 1)], field[i + 1 + ni * j], field[i + 1 + ni * (j + 1)], fx, fy);
}

glm::vec2 InterpolateGradientOnGrid(const glm::vec2& pos_normalized, const std::vector<float>& field, const int ni, const int nj) {
	int i, j;
	float fx, fy;

	GetBaryCentric(pos_normalized.x, i, fx, 0, ni - 2);
	GetBaryCentric(pos_normalized.y, j, fy, 0, nj - 2);
	
	float v00 = field[i + ni * j];
	float v01 = field[i + (ni + 1) * j];
	float v10 = field[i + 1 + ni * j];
	float v11 = field[i + 1 + (ni + 1) * j];

	float dx0 = v10 - v00;
	float dx1 = v11 - v01;

	float dv_dx = Lerp(dx0, dx1, fy);

	float dy0 = v01 - v00;
	float dy1 = v11 - v10;

	float dv_dy = Lerp(dy0, dy1, fx);

	return glm::vec2(dv_dx, dv_dy);
}

float FractionInside(float sdf_left, float sdf_right) {
	if (sdf_left < 0.0f && sdf_right < 0.0f) {
		return 1.0f;
	}
	if (sdf_left > 0.0f && sdf_right > 0.0f) {
		return 0.0f;
	}
	if (sdf_left > 0.0f && sdf_right < 0.0f) {
		return sdf_right / sdf_right - sdf_left;
	}
	if (sdf_left < 0.0f && sdf_right > 0.0f) {
		return sdf_left / sdf_left - sdf_right;
	}
	return 0.0f;
}
