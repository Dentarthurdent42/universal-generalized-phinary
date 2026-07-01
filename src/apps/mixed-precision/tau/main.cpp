/** ************************************************************************
* tau: demonstration of the balanced tau (phi^2) number system
*
* The balanced tau number system is the "ternary phi-representation" of the
* Generalized Phi Number System (GPNS) by S. Agaian and Y. Zhou,
* "Generalized Phi Number System and its Applications for Image Decomposition
* and Enhancement", Proc. SPIE-IS&T Electronic Imaging, Vol. 7881, 2011.
*
* It is a positional number system with the irrational base
*     tau = phi^2 = (3 + sqrt(5)) / 2  ~= 2.6180339887...
* and the balanced-ternary digit set { -1, 0, +1 } (printed here as -, 0, +).
*
* @copyright:  Copyright (c) 2017 Stillwater Supercomputing, Inc.
* @license:    MIT Open Source license
*
* This file is part of the MPADAO project
* *************************************************************************
*/
#include <iostream>
#include <iomanip>

#include <universal/number/tau/tau.hpp>

int main()
try {
	using namespace sw::universal;

	std::cout << "Balanced tau (phi^2) number system\n";
	std::cout << "tau = phi^2 = (3 + sqrt(5))/2 ~= 2.6180339887...\n\n";

	using Tau = tau<16, 6>;   // 10 integer trits, 6 fractional trits

	std::cout << "configuration : " << type_tag<Tau>() << '\n';
	std::cout << "fields        : " << type_field<Tau>() << '\n';
	std::cout << "range         : " << range<Tau>() << "\n\n";

	// the defining algebraic identities of base tau
	std::cout << "Defining identities (tau^2 = 3*tau - 1, so tau + 1/tau = 3):\n";
	{
		Tau t; t = Tau::tau_base;
		std::cout << "  tau         = " << std::setw(10) << t.to_double() << "   " << to_triplets(t) << '\n';
		Tau three = t + Tau(1.0 / Tau::tau_base);
		std::cout << "  tau + 1/tau = " << std::setw(10) << three.to_double() << "   " << to_triplets(three) << '\n';
	}
	std::cout << '\n';

	// the GPNS ternary phi-representation of the first integers (paper Table 1)
	std::cout << "Balanced-tau representation of integers:\n";
	std::cout << "   n      value   trits (- = -1, + = +1)\n";
	for (int n = 0; n <= 12; ++n) {
		Tau v = n;
		std::cout << std::setw(4) << n
			<< std::setw(11) << v.to_double()
			<< "   " << to_triplets(v) << '\n';
	}
	std::cout << '\n';

	// arithmetic
	std::cout << "Arithmetic (rounded to the nearest representable tau value):\n";
	Tau a = 2.0, b = 3.0;
	std::cout << "  a = " << a << ", b = " << b << '\n';
	std::cout << "  a + b = " << (a + b) << "   " << to_triplets(a + b) << '\n';
	std::cout << "  a - b = " << (a - b) << "   " << to_triplets(a - b) << '\n';
	std::cout << "  a * b = " << (a * b) << "   " << to_triplets(a * b) << '\n';
	std::cout << "  b / a = " << (b / a) << "   " << to_triplets(b / a) << '\n';
	std::cout << '\n';

	// symmetric range is a hallmark of a balanced digit set
	std::cout << "Symmetry: maxneg == -maxpos ? "
		<< std::boolalpha
		<< (Tau(SpecificValue::maxneg) == -Tau(SpecificValue::maxpos)) << '\n';

	return EXIT_SUCCESS;
}
catch (const std::exception& e) {
	std::cerr << "Caught exception: " << e.what() << '\n';
	return EXIT_FAILURE;
}
