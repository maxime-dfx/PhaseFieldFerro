#pragma once
#include <string>
#include <vector>
#include <functional>
#include <iostream>
#include <cmath>
#include <exception>
#include <sstream>

// Testing
// -----------------------------------------------------------------------
// Framework de test minimal, header-only, sans dependance externe (pas de
// GoogleTest/Catch2 - le sandbox de build n'a pas d'acces reseau garanti
// pour les recuperer). Suffisant pour des tests unitaires/regression avec
// enregistrement automatique via macro, un peu comme GTest.
//
// Usage :
//   TEST_CASE("mon_test") {
//       CHECK_NEAR(2.0 + 2.0, 4.0, 1e-12);
//       CHECK_TRUE(1 < 2);
//   }
//
// Toutes les TEST_CASE de tous les .cpp lies dans l'executable de test
// s'enregistrent automatiquement (construction statique) ; il suffit
// d'appeler Testing::run_all() depuis un main().
namespace Testing {

class TestFailure : public std::exception {
public:
    std::string msg;
    explicit TestFailure(std::string m) : msg(std::move(m)) {}
    const char* what() const noexcept override { return msg.c_str(); }
};

struct TestCaseEntry {
    std::string name;
    std::function<void()> fn;
};

inline std::vector<TestCaseEntry>& registry() {
    static std::vector<TestCaseEntry> tests;
    return tests;
}

struct Registrar {
    Registrar(const std::string& name, std::function<void()> fn) {
        registry().push_back({name, std::move(fn)});
    }
};

// Renvoie 0 si tout passe, 1 sinon (utilisable directement comme code de
// retour de main() - conventionnel pour l'integration CTest/CI).
inline int run_all() {
    int passed = 0, failed = 0;
    std::cout << "=== PhaseFieldFerro : " << registry().size() << " tests ===\n\n";
    for (auto& t : registry()) {
        try {
            t.fn();
            std::cout << "[PASS] " << t.name << "\n";
            ++passed;
        } catch (const TestFailure& e) {
            std::cout << "[FAIL] " << t.name << "\n       " << e.msg << "\n";
            ++failed;
        } catch (const std::exception& e) {
            std::cout << "[FAIL] " << t.name << " (exception non prevue)\n       " << e.what() << "\n";
            ++failed;
        }
    }
    std::cout << "\n" << passed << " reussis, " << failed << " echoues, sur " << registry().size() << " tests.\n";
    return failed == 0 ? 0 : 1;
}

} // namespace Testing

#define TEST_CONCAT_(a,b) a##b
#define TEST_CONCAT(a,b) TEST_CONCAT_(a,b)

// Declare un cas de test. Le nom (string) doit etre unique et descriptif -
// il sert d'identifiant dans le rapport, pas le nom de la fonction C++.
#define TEST_CASE(name) \
    static void TEST_CONCAT(test_fn_, __LINE__)(); \
    static ::Testing::Registrar TEST_CONCAT(test_reg_, __LINE__)(name, TEST_CONCAT(test_fn_, __LINE__)); \
    static void TEST_CONCAT(test_fn_, __LINE__)()

#define CHECK_NEAR(actual, expected, tol) \
    do { \
        double _fdd_a = static_cast<double>(actual); \
        double _fdd_e = static_cast<double>(expected); \
        double _fdd_t = static_cast<double>(tol); \
        if (!std::isfinite(_fdd_a) || std::abs(_fdd_a - _fdd_e) > _fdd_t) { \
            std::ostringstream _fdd_oss; \
            _fdd_oss << #actual " != " #expected " : " << _fdd_a << " vs " << _fdd_e \
                     << " (tol=" << _fdd_t << ") a " __FILE__ ":" << __LINE__; \
            throw ::Testing::TestFailure(_fdd_oss.str()); \
        } \
    } while (0)

#define CHECK_TRUE(cond) \
    do { \
        if (!(cond)) { \
            throw ::Testing::TestFailure(std::string("Condition fausse : " #cond " a " __FILE__ ":") + std::to_string(__LINE__)); \
        } \
    } while (0)
