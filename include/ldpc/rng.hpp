/**
 * @file rng.hpp
 * @brief Defines the ldpc::rng::RandomNumberGenerator class, which generates random double values between 0 and 1
 */

#ifndef RNG_HPP
#define RNG_HPP

#include <random>
#include <chrono> // for std::chrono::system_clock
#include <algorithm>

namespace ldpc {
    namespace rng {

        /**
         * @brief Generates random double values between 0 and 1 using a Mersenne Twister random number generator
         *
         * The RandomNumberGenerator class creates a Mersenne Twister random number generator, and provides a method
         * for generating random double values between 0 and 1 using a uniform distribution.
         */
        class RandomNumberGenerator {
        public:
            /**
             * @brief Constructs a new RandomNumberGenerator object with an optional seed
             *
             * If no seed is specified, the generator is seeded using the system clock.
             *
             * @param seed The seed value to use for the random number generator
             */
            explicit RandomNumberGenerator(int seed = 0) : gen(seed) {
                if (seed == 0) {
                    gen.seed(std::chrono::system_clock::now().time_since_epoch().count());
                }
            }

            ~RandomNumberGenerator() = default;

            /**
             * @brief Generates a new random double value between 0 and 1
             *
             * @return A random double value between 0 and 1
             */
            double random_double() {
                // Create a uniform distribution between 0 and 1
                std::uniform_real_distribution<double> dis(0.0, 1.0);

                // Generate a random number between 0 and 1
                double random_num = dis(gen);

                return random_num;
            }

            /**
             * @brief Generates a new random integer value between 0 and max_int
             *
             * @param max_int The maximum value that the random integer can take
             * @return A random integer value between 0 and max_int
             */
            int random_int(int max_int) {
                // Create a uniform distribution between 0 and max_int
                std::uniform_int_distribution<int> dis(0, max_int);

                // Generate a random number between 0 and max_int
                int random_num = dis(gen);

                return random_num;
            }

        private:
            std::mt19937 gen; /**< The Mersenne Twister random number generator used by the class */
        };


        /**
         * @brief A templated class for shuffling lists of data.
         *
         * This class provides the capability to shuffle lists of data of any type using a
         * random number generator.
         *
         * @tparam T The type of data to shuffle.
         */
        template<typename T>
        class RandomListShuffle {
        public:
            /**
             * @brief Default constructor.
             */
            RandomListShuffle() = default;

            /**
             * @brief Constructor that allows specifying a seed.
             *
             * This constructor initializes the random number generator with the provided seed.
             * If the seed is zero, it uses the system clock to generate a seed.
             *
             * @param seed The seed for the random number generator.
             */
            explicit RandomListShuffle(unsigned int seed) {
                this->seed(seed);
            }

            /**
             * @brief Set the seed for the random number generator.
             *
             * This function allows you to set a custom seed for the random number generator.
             * If the seed is zero, it uses the system clock to generate a seed.
             *
             * @param seed The seed for the random number generator.
             */
            void seed(unsigned int seed) {
                if (seed == 0) {
                    // Use the system clock to generate a seed
                    auto now = std::chrono::system_clock::now();
                    seed = static_cast<unsigned int>(now.time_since_epoch().count());
                }
                generator.seed(seed);
            }

            /**
             * @brief Shuffle a vector of data.
             *
             * This function shuffles the elements in the provided vector using the random
             * number generator.
             *
             * @param data The vector of data to be shuffled.
             */
            void shuffle(std::vector<T> &data) {
                std::shuffle(data.begin(), data.end(), generator);
            }

        private:
            std::mt19937 generator; ///< The random number generator.
        };

    }
}  // namespace ldpc::rng

#endif // RNG_HPP
