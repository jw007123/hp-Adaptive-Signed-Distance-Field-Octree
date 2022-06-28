#pragma once

#include "Literals.h"

/*
	Functions and Lookups
*/
namespace SDF
{
	constexpr usize BASIS_MAX_DEGREE = 12;
	constexpr usize TREE_MAX_DEPTH   = 10;

	/*
		Constexpr-satisfying versions of Pow and Sqrt.
	*/
	constexpr f64 PowConst(const f64 base_, const usize pow_)
	{
		if (pow_ == 0)
		{
			return 1.0;
		}
		else
		{
			return base_ * PowConst(base_, pow_ - 1);
		}
	}
	constexpr f64 SqrtConst(const f64 x_)
	{
		// Newton's Method
		f64 guess = x_;
		for (usize i = 0; i < 100; ++i)
		{
			guess = 0.5 * (guess + x_ / guess);
		}

		return guess;
	}

	/*
		Mod3Lookup[i][j] = (i + j) % 3. Used to determine which bits of a node to use when integrating over shared faces.
	*/
	u8 Mod3Lookup[3][3] =
	{
		{0, 1, 2},
		{1, 2, 0},
		{2, 0, 1}
	};

	/*
		Standard sum to N. store[n] = n * (n + 1) / 2
	*/
	struct SumToNCalc
	{
		constexpr SumToNCalc() : values()
		{
			for (u32 i = 0; i <= (4 * BASIS_MAX_DEGREE + 1); ++i)
			{
                u32 sum = 0;
				for (u32 j = 0; j <= i; ++j)
				{
					sum += j;
				}

				values[i] = sum;
			}
		}
        u32 values[(4 * BASIS_MAX_DEGREE + 1) + 1];
	};
	constexpr const u32(&SumToN)[(4 * BASIS_MAX_DEGREE + 1) + 1] = SumToNCalc().values;

	/*
		In equation (4), precomputes the sqrt portions under the assumption that the internal octree
		volume is taken to be the unit cube.
	*/
	struct NormalisedLengthsCalc
	{
		constexpr NormalisedLengthsCalc() : values()
		{
			for (u32 i = 0; i <= BASIS_MAX_DEGREE; ++i)
			{
				for (u32 j = 0; j <= TREE_MAX_DEPTH; ++j)
				{
					// Sizes are 1 at level 0, 0.5 at level 1, ... , 2^-n at level n
					values[i][j] = SqrtConst((2.0 * i + 1.0) * PowConst(2.0, j));
				}
			}
		}
		f64 values[BASIS_MAX_DEGREE + 1][TREE_MAX_DEPTH + 1];
	};
	constexpr const f64(&NormalisedLengths)[BASIS_MAX_DEGREE + 1][TREE_MAX_DEPTH + 1] = NormalisedLengthsCalc().values;

	/*
		Stores the number of coefficient in each basis polynomial. I.e. store[n] tells us that
		the nth Legendre polynomial contains store[n] coefficients. If the dimension is M and
		the max basis degree is N, then the formula is:
				
					   store_M^N[n] = 1 / M! * (N + 1) * (N + 2) * ... * (N + M)
	*/
	struct LegendreCoefficientCountCalc
	{
		constexpr LegendreCoefficientCountCalc() : values()
		{
			const f64 f = 1.0 / 6.0;

			for (u32 i = 0; i <= BASIS_MAX_DEGREE; ++i)
			{
				values[i] = (u32)(f * (i + 1) * (i + 2) * (i + 3));
			}
		}
        u32 values[BASIS_MAX_DEGREE + 1];

		static constexpr const u32 Size()
		{
			const f64 f = 1.0 / 6.0;
			return (u32)(f * (BASIS_MAX_DEGREE + 1) * (BASIS_MAX_DEGREE + 2) * (BASIS_MAX_DEGREE + 3));
		}
	};
	constexpr const u32(&LegendreCoeffientCount)[BASIS_MAX_DEGREE + 1] = LegendreCoefficientCountCalc().values;

	/*
		Stores the constants used in each term of the reccurence relation definition for Legendre polynomials.
		Removes the need to perform repeated divides when querying the tree.
	*/
	struct LegendreCoefficientCalc
	{
		constexpr LegendreCoefficientCalc() : values()
		{
			values[0][0] = 0.0;
			values[0][1] = 0.0;

			for (u32 i = 1; i <= BASIS_MAX_DEGREE; ++i)
			{
				values[i][0] = (2.0 * i - 1.0) / i;
				values[i][1] = (i - 1.0) / i;
			}
		}
		f64 values[BASIS_MAX_DEGREE + 1][2];
	};
	constexpr const f64(&LegendreCoefficent)[BASIS_MAX_DEGREE + 1][2] = LegendreCoefficientCalc().values;

	/*
		Stores all possible combinations of BasisIndex for dimension M and max basis degree N. Saves having
		to compute these on the fly or store them inside the node (both very bad ideas!).
	*/
	struct BasisIndexValuesCalc
	{
		constexpr BasisIndexValuesCalc() : values()
		{
            u32 valuesIdx = 0;
			for (u32 p = 0; p <= BASIS_MAX_DEGREE; ++p)
			{
				for (u32 i = 0; i <= p; ++i)
				{
					for (u32 j = 0; j <= p - i; ++j)
					{
						for (u32 k = 0; k <= p - i - j; ++k)
						{
							if ((i + j + k) == p)
							{
								values[valuesIdx][0] = i;
								values[valuesIdx][1] = j;
								values[valuesIdx][2] = k;
								valuesIdx++;
							}
						}
					}
				}
			}
		}
        u32 values[LegendreCoefficientCountCalc::Size()][3];
	};
	constexpr const u32(&BasisIndexValues)[LegendreCoefficientCountCalc::Size()][3] = BasisIndexValuesCalc().values;
}