#include "Definitions.hpp"
#include "QuantumComputation.hpp"
#include "dd/Package.hpp"
#include "dd/Simulation.hpp"
#include "ecc/Id.hpp"
#include "ecc/Q18Surface.hpp"
#include "ecc/Q3Shor.hpp"
#include "ecc/Q5Laflamme.hpp"
#include "ecc/Q7Steane.hpp"
#include "ecc/Q9Shor.hpp"
#include "ecc/Q9Surface.hpp"
#include "operations/NonUnitaryOperation.hpp"
#include "operations/OpType.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <gtest/gtest.h>
#include <iostream>
#include <memory>
#include <vector>

using namespace qc;

using circuitFunctions =
        std::function<std::shared_ptr<qc::QuantumComputation>()>;

struct TestCase {
    circuitFunctions circuit;
    bool             testNoise;
    size_t           insertNoiseAfterNGates;
};

class DDECCFunctionalityTest : public ::testing::Test {
protected:
    void SetUp() override {}

    void TearDown() override {}

    static std::shared_ptr<qc::QuantumComputation> createIdentityCircuit() {
        auto qc = std::make_shared<qc::QuantumComputation>();
        qc->addQubitRegister(1U);
        qc->addClassicalRegister(1U, "resultReg");
        qc->x(0);
        qc->x(0);
        qc->measure(0, {"resultReg", 0});
        return qc;
    }

    static std::shared_ptr<qc::QuantumComputation> createXCircuit() {
        auto qc = std::make_shared<qc::QuantumComputation>();
        qc->addQubitRegister(1U);
        qc->addClassicalRegister(1U, "resultReg");
        qc->x(0);
        qc->measure(0, {"resultReg", 0});
        return qc;
    }

    static std::shared_ptr<qc::QuantumComputation> createYCircuit() {
        auto qc = std::make_shared<qc::QuantumComputation>();
        qc->addQubitRegister(1U);
        qc->addClassicalRegister(1U, "resultReg");
        qc->h(0);
        qc->y(0);
        qc->h(0);
        qc->measure(0, {"resultReg", 0});
        return qc;
    }

    static std::shared_ptr<qc::QuantumComputation> createHCircuit() {
        auto qc = std::make_shared<qc::QuantumComputation>();
        qc->addQubitRegister(1U);
        qc->addClassicalRegister(1U, "resultReg");
        qc->h(0);
        qc->measure(0, {"resultReg", 0});
        return qc;
    }

    static std::shared_ptr<qc::QuantumComputation> createHTCircuit() {
        auto qc = std::make_shared<qc::QuantumComputation>();
        qc->addQubitRegister(1U);
        qc->addClassicalRegister(1U, "resultReg");
        qc->h(0);
        qc->t(0);
        qc->tdg(0);
        qc->h(0);
        qc->measure(0, {"resultReg", 0});
        return qc;
    }

    static std::shared_ptr<qc::QuantumComputation> createXSCircuit() {
        auto qc = std::make_shared<qc::QuantumComputation>();
        qc->addQubitRegister(1U);
        qc->addClassicalRegister(1U, "resultReg");
        qc->x(0);
        qc->s(0);
        qc->sdg(0);
        qc->x(0);
        qc->measure(0, {"resultReg", 0});
        return qc;
    }

    static std::shared_ptr<qc::QuantumComputation> createHZCircuit() {
        auto qc = std::make_shared<qc::QuantumComputation>();
        qc->addQubitRegister(1U);
        qc->addClassicalRegister(1U, "resultReg");
        qc->h(0);
        qc->z(0);
        qc->h(0);
        qc->measure(0, {"resultReg", 0});
        return qc;
    }

    static std::shared_ptr<qc::QuantumComputation> createCXCircuit() {
        auto qc = std::make_shared<qc::QuantumComputation>();
        qc->addQubitRegister(2U);
        qc->addClassicalRegister(2U, "resultReg");
        qc->x(0);
        qc->cx(0, 1);
        qc->measure(0, {"resultReg", 0});
        qc->measure(1, {"resultReg", 1});
        return qc;
    }

    static std::shared_ptr<qc::QuantumComputation> createCZCircuit() {
        auto qc = std::make_shared<qc::QuantumComputation>();
        qc->addQubitRegister(2U);
        qc->addClassicalRegister(2U, "resultReg");
        qc->x(0);
        qc->h(1);
        qc->cz(0, 1);
        qc->h(1);
        qc->measure(0, {"resultReg", 0});
        qc->measure(1, {"resultReg", 1});
        return qc;
    }

    static std::shared_ptr<qc::QuantumComputation> createCYCircuit() {
        auto qc = std::make_shared<qc::QuantumComputation>();
        qc->addQubitRegister(2U);
        qc->addClassicalRegister(2U, "resultReg");
        qc->x(0);
        qc->cy(0, 1);
        qc->measure(0, {"resultReg", 0});
        qc->measure(1, {"resultReg", 1});
        return qc;
    }

    template <class eccType>
    [[nodiscard]] bool
    testCircuits(const std::vector<TestCase>& circuitsExpectToPass) const {
        size_t circuitCounter = 0;
        for (const auto& testParameter : circuitsExpectToPass) {
            auto qcOriginal = testParameter.circuit();
            auto mapper     = std::make_unique<eccType>(qcOriginal, 0);
            mapper->apply();
            circuitCounter++;
            std::cout << "Testing circuit " << circuitCounter << "\n";
            bool const success = testErrorCorrectionCircuit(
                    mapper->getOriginalCircuit(), mapper->getMappedCircuit(),
                    testParameter.testNoise, mapper->getDataQubits(),
                    testParameter.insertNoiseAfterNGates);
            if (!success) {
                return false;
            }
        }
        return true;
    }

    static bool testErrorCorrectionCircuit(
            const std::shared_ptr<qc::QuantumComputation>& qcOriginal,
            const std::shared_ptr<qc::QuantumComputation>& qcMapped,
            bool simulateWithErrors, const std::vector<Qubit>& dataQubits = {},
            const std::size_t insertErrorAfterNGates = 0) {
        if (!simulateWithErrors) {
            return simulateAndVerifyResults(qcOriginal, qcMapped);
        }
        return (std::all_of(
                dataQubits.begin(), dataQubits.end(),
                [qcMapped, insertErrorAfterNGates, qcOriginal,
                 simulateWithErrors](Qubit target) {
                    auto initialQuit = qcMapped->begin();
                    qcMapped->insert(
                            initialQuit + static_cast<int64_t>(insertErrorAfterNGates),
                            std::make_unique<NonUnitaryOperation>(std::vector<Qubit>{target},
                                                                  qc::Reset));
                    auto result =
                            simulateAndVerifyResults(qcOriginal, qcMapped, simulateWithErrors,
                                                     insertErrorAfterNGates, target);
                    qcMapped->erase(initialQuit +
                                    static_cast<int64_t>(insertErrorAfterNGates));
                    return result;
                }));
    }

    static bool simulateAndVerifyResults(
            const std::shared_ptr<qc::QuantumComputation>& qcOriginal,
            const std::shared_ptr<qc::QuantumComputation>& qcMapped,
            bool simulateWithErrors = false, std::size_t insertErrorAfterNGates = 0,
            Qubit target = 0, double tolerance = 0.2, std::size_t shots = 50,
            std::size_t seed = 5) {
        auto toleranceAbsolute =
                (static_cast<double>(shots) / 100.0) * (tolerance * 100.0);

        auto ddOriginal       = std::make_unique<dd::Package<>>(qcOriginal->getNqubits());
        auto originalRootEdge = ddOriginal->makeZeroState(qcOriginal->getNqubits());
        ddOriginal->incRef(originalRootEdge);

        auto measurementsOriginal =
                simulate(qcOriginal.get(), originalRootEdge, *ddOriginal, shots);

        auto ddEcc       = std::make_unique<dd::Package<>>(qcMapped->getNqubits());
        auto eccRootEdge = ddEcc->makeZeroState(qcMapped->getNqubits());
        ddEcc->incRef(eccRootEdge);

        auto measurementsProtected =
                simulate(qcMapped.get(), eccRootEdge, *ddEcc, shots, seed);
        for (auto const& [classicalBit, hits] : measurementsOriginal) {
            // Since the result is stored as one bit string. I have to count the
            // relevant classical bits.
            size_t eccHits = 0;
            for (auto const& [eccMeasure, tempHits] : measurementsProtected) {
                if (0 == eccMeasure.compare(eccMeasure.length() - classicalBit.length(),
                                            classicalBit.length(), classicalBit)) {
                    eccHits += tempHits;
                }
            }
            auto difference = std::max(eccHits, hits) - std::min(eccHits, hits);
            std::cout << "Diff/tolerance " << difference << "/" << toleranceAbsolute
                      << " Original register: " << hits
                      << " ecc register: " << eccHits << "\n";
            if (simulateWithErrors) {
                std::cout << " Simulating an error in qubit " << target << " after "
                          << insertErrorAfterNGates << " gates.\n";
            }
            if (static_cast<double>(difference) > toleranceAbsolute) {
                std::cout << "Error is too large!\n";
                return false;
            }
        }
        return true;
    }
};

TEST_F(DDECCFunctionalityTest, testIdEcc) {
    const size_t insertNoiseAfterNQubits = 1;

    std::vector<TestCase> const circuitsExpectToPass = {
            {createIdentityCircuit, false, insertNoiseAfterNQubits},
            {createXCircuit, false, insertNoiseAfterNQubits},
            {createYCircuit, false, insertNoiseAfterNQubits},
            {createHCircuit, false, insertNoiseAfterNQubits},
            {createHTCircuit, false, insertNoiseAfterNQubits},
            {createXSCircuit, false, insertNoiseAfterNQubits},
            {createHZCircuit, false, insertNoiseAfterNQubits},
            {createCXCircuit, false, insertNoiseAfterNQubits},
            {createCZCircuit, false, insertNoiseAfterNQubits},
            {createCYCircuit, false, insertNoiseAfterNQubits},
    };
    EXPECT_TRUE(testCircuits<ecc::Id>(circuitsExpectToPass));
}

TEST_F(DDECCFunctionalityTest, testQ3Shor) {
    const size_t insertNoiseAfterNQubits = 4;

    std::vector<TestCase> const circuitsExpectToPass = {
            {createIdentityCircuit, true, insertNoiseAfterNQubits},
            {createXCircuit, true, insertNoiseAfterNQubits},
            {createXSCircuit, true, insertNoiseAfterNQubits},
            {createCXCircuit, true, insertNoiseAfterNQubits * 2},
            {createCYCircuit, true, insertNoiseAfterNQubits * 2},
    };
    std::vector<TestCase> const circuitsExpectToFail = {
            {createYCircuit, true, insertNoiseAfterNQubits},
            {createHCircuit, true, insertNoiseAfterNQubits},
            {createHTCircuit, true, insertNoiseAfterNQubits},
            {createCZCircuit, true, insertNoiseAfterNQubits},
            {createHZCircuit, true, insertNoiseAfterNQubits},
    };
    EXPECT_TRUE(testCircuits<ecc::Q3Shor>(circuitsExpectToPass));
    EXPECT_ANY_THROW([[maybe_unused]] const bool result =
                             testCircuits<ecc::Q3Shor>(circuitsExpectToFail));
}

TEST_F(DDECCFunctionalityTest, testQ5LaflammeEcc) {
    const size_t insertNoiseAfterNQubits = 61;

    std::vector<TestCase> const circuitsExpectToPass = {
            {createIdentityCircuit, true, insertNoiseAfterNQubits},
            {createXCircuit, true, insertNoiseAfterNQubits},
    };
    std::vector<TestCase> const circuitsExpectToFail = {
            {createYCircuit, true, insertNoiseAfterNQubits},
            {createHCircuit, true, insertNoiseAfterNQubits},
            {createHTCircuit, true, insertNoiseAfterNQubits},
            {createHZCircuit, true, insertNoiseAfterNQubits},
            {createXSCircuit, true, insertNoiseAfterNQubits},
            {createCXCircuit, true, insertNoiseAfterNQubits},
            {createCZCircuit, true, insertNoiseAfterNQubits},
            {createCYCircuit, true, insertNoiseAfterNQubits},
    };
    EXPECT_TRUE(testCircuits<ecc::Q5Laflamme>(circuitsExpectToPass));
    EXPECT_ANY_THROW([[maybe_unused]] const bool result =
                             testCircuits<ecc::Q5Laflamme>(circuitsExpectToFail));
}

TEST_F(DDECCFunctionalityTest, testQ7Steane) {
    const size_t insertNoiseAfterNQubits = 57;

    std::vector<TestCase> const circuitsExpectToPass = {
            {createIdentityCircuit, true, insertNoiseAfterNQubits},
            {createXCircuit, true, insertNoiseAfterNQubits},
            {createYCircuit, true, insertNoiseAfterNQubits},
            {createHCircuit, true, insertNoiseAfterNQubits},
            {createHTCircuit, true, insertNoiseAfterNQubits},
            {createHZCircuit, true, insertNoiseAfterNQubits},
            {createXSCircuit, true, insertNoiseAfterNQubits},
            {createCXCircuit, true, insertNoiseAfterNQubits * 2},
            {createCZCircuit, true, insertNoiseAfterNQubits * 2},
            {createCYCircuit, true, insertNoiseAfterNQubits * 2},
    };
    EXPECT_TRUE(testCircuits<ecc::Q7Steane>(circuitsExpectToPass));
}

TEST_F(DDECCFunctionalityTest, testQ9ShorEcc) {
    const size_t insertNoiseAfterNQubits = 7;

    std::vector<TestCase> const circuitsExpectToPass = {
            {createIdentityCircuit, true, insertNoiseAfterNQubits},
            {createXCircuit, true, insertNoiseAfterNQubits},
            {createCXCircuit, true, insertNoiseAfterNQubits * 2},
            {createCYCircuit, true, insertNoiseAfterNQubits * 2},
    };

    std::vector<TestCase> const circuitsExpectToFail = {
            {createYCircuit, true, insertNoiseAfterNQubits},
            {createHCircuit, true, insertNoiseAfterNQubits},
            {createHTCircuit, true, insertNoiseAfterNQubits},
            {createHZCircuit, true, insertNoiseAfterNQubits},
            {createXSCircuit, true, insertNoiseAfterNQubits},
            {createCZCircuit, true, insertNoiseAfterNQubits},
    };

    EXPECT_TRUE(testCircuits<ecc::Q9Shor>(circuitsExpectToPass));
    EXPECT_ANY_THROW([[maybe_unused]] const bool result =
                             testCircuits<ecc::Q9Shor>(circuitsExpectToFail));
}

TEST_F(DDECCFunctionalityTest, testQ9SurfaceEcc) {
    const size_t insertNoiseAfterNQubits = 55;

    std::vector<TestCase> const circuitsExpectToPass = {
            {createIdentityCircuit, true, insertNoiseAfterNQubits},
            {createXCircuit, true, insertNoiseAfterNQubits},
            {createYCircuit, true, insertNoiseAfterNQubits},
            {createHCircuit, true, insertNoiseAfterNQubits},
            {createHZCircuit, true, insertNoiseAfterNQubits},
    };

    std::vector<TestCase> const circuitsExpectToFail = {
            {createHTCircuit, true, insertNoiseAfterNQubits},
            {createXSCircuit, true, insertNoiseAfterNQubits},
            {createCXCircuit, true, insertNoiseAfterNQubits},
            {createCZCircuit, true, insertNoiseAfterNQubits},
            {createCYCircuit, true, insertNoiseAfterNQubits},
    };

    EXPECT_TRUE(testCircuits<ecc::Q9Surface>(circuitsExpectToPass));

    EXPECT_ANY_THROW([[maybe_unused]] const bool result =
                             testCircuits<ecc::Q9Surface>(circuitsExpectToFail));
}

TEST_F(DDECCFunctionalityTest, testQ18SurfaceEcc) {
    const size_t insertNoiseAfterNQubits = 47;

    std::vector<TestCase> const circuitsExpectToPass{
            {createIdentityCircuit, false, insertNoiseAfterNQubits},
            {createXCircuit, false, insertNoiseAfterNQubits},
            {createYCircuit, false, insertNoiseAfterNQubits},
            {createHCircuit, false, insertNoiseAfterNQubits},
            {createHZCircuit, false, insertNoiseAfterNQubits},
    };

    std::vector<TestCase> const circuitsExpectToFail = {
            {createHTCircuit, true, insertNoiseAfterNQubits},
            {createXSCircuit, true, insertNoiseAfterNQubits},
            {createCXCircuit, true, insertNoiseAfterNQubits},
            {createCZCircuit, true, insertNoiseAfterNQubits},
            {createCYCircuit, true, insertNoiseAfterNQubits},
    };

    EXPECT_TRUE(testCircuits<ecc::Q18Surface>(circuitsExpectToPass));
    EXPECT_ANY_THROW([[maybe_unused]] const bool result =
                             testCircuits<ecc::Q18Surface>(circuitsExpectToFail));
}
