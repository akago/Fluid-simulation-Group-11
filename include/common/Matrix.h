#pragma once

#include <cmath>
#include <gfx/vec2.h>

struct matrix {
	Vec2f row1;
	Vec2f row2;
	
	matrix() : row1(Vec2f(0.0, 0.0)), row2(Vec2f(0.0, 0.0)) {}
	matrix(Vec2f r1, Vec2f r2) : row1(r1), row2(r2) {}

	matrix multiByMat(matrix other) {
		matrix result;
		result.row1[0] = row1[0] * other.row1[0] + row1[1] * other.row2[0];
		result.row1[1] = row1[0] * other.row1[1] + row1[1] * other.row2[1];
		result.row2[0] = row2[0] * other.row1[0] + row2[1] * other.row2[0];
		result.row2[1] = row2[0] * other.row1[1] + row2[1] * other.row2[1];
		return result;
	}

	Vec2f multiByVec2(Vec2f vec) {
		Vec2f result;

		result[0] = row1 * vec;
		result[1] = row2 * vec;
		
		return result;
	}
};
