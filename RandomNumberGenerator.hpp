#pragma once
#include <random>
#include <limits>
#include <type_traits>
#include <algorithm>
#include <cassert>
#include <iostream>

/*
EXAMPLE
		RandomNumberGenerator<>::instance().seed(1);
		double r[10]; // 10 doubles randomicos entre 0 e 1
		for (int i = 0; i < 10; ++i) {
			r[i] = random_number<double>(0.0, 1.0);
			// se for int, poderia ser por exemplo
			// r[i] = random_number<int>(0, 25000)
			cout << r[i] << " ";
		}
		cout << endl;
		getchar();
 */

template<typename RandomEngine = std::mt19937>
class RandomNumberGenerator {
public:
    RandomNumberGenerator() : engine_(device_()) { }

	RandomNumberGenerator(typename RandomEngine::result_type initialSeed) 
		: engine_(initialSeed) {
	}

    void randomize() {
        engine_.seed(device_());
    }

	void seed(typename RandomEngine::result_type value) {
		engine_.seed(value);
	}

    template<typename T>
    typename std::enable_if<std::is_integral<T>::value, T>::type
        generate(T from, T to) {
        static std::uniform_int_distribution<T> d;
        using parm_t = typename decltype(d)::param_type;
        return d(engine_, parm_t {from, to});
    }

    template<typename T>
    typename std::enable_if<std::is_floating_point<T>::value, T>::type
        generate(T from, T to) {
        static std::uniform_real_distribution<T> d;
        using parm_t = typename decltype(d)::param_type;
        return d(engine_, parm_t {from, to});
    }

    static RandomNumberGenerator& instance() {
        static RandomNumberGenerator rng;
        return rng;
    }

	RandomEngine& engine() { return engine_; }
	std::random_device& device() { return device_; }

private:
    std::random_device device_;
    RandomEngine engine_;
};

template<typename T>
inline T random_number(T from = std::numeric_limits<T>::min(),
                      T to = std::numeric_limits<T>::max()) {
    return RandomNumberGenerator<>::instance().generate(from, to);
}

template<typename StringType>
inline StringType random_string(int size, const StringType& range) {    
    int rangeSize = range.size();
    StringType s;
    if (rangeSize > 0) {
        s.resize(size);
        std::generate_n(std::begin(s), size, [&](){ 
            return range[random_number(0, rangeSize - 1)];
        });
    }     
    return s;
}
