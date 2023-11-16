#pragma once
#include <vector>

class PPEMatrix2D
{
public:
	unsigned int n_;
	std::vector<double> center_;
	std::vector<double> up_;
	std::vector<double> right_;

	explicit PPEMatrix2D(unsigned int size);

	void resize(int size, int ni, int nj) {
		center_.resize(size);
		up_.resize(size - )
	}
private:


};
