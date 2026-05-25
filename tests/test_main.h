#pragma once

#include <Eigen/Dense>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace test {

inline void fail(const char* expr, const char* file, int line, const std::string& message = "")
{
    std::ostringstream out;
    out << file << ":" << line << ": requirement failed: " << expr;
    if (!message.empty()) out << " (" << message << ")";
    throw std::runtime_error(out.str());
}

inline bool near(double a, double b, double tol)
{
    return std::abs(a - b) <= tol;
}

inline bool nearVec(const Eigen::VectorXd& a, const Eigen::VectorXd& b, double tol)
{
    return (a - b).norm() <= tol;
}

inline bool finiteVec(const Eigen::VectorXd& v)
{
    return v.array().isFinite().all();
}

template <typename Fn>
int run(const char* name, Fn&& fn)
{
    try {
        fn();
        std::cout << "[PASS] " << name << "\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& e) {
        std::cerr << "[FAIL] " << name << ": " << e.what() << "\n";
        return EXIT_FAILURE;
    }
}

} // namespace test

#define REQUIRE(expr) do { if (!(expr)) ::test::fail(#expr, __FILE__, __LINE__); } while (false)
#define REQUIRE_MSG(expr, msg) do { if (!(expr)) ::test::fail(#expr, __FILE__, __LINE__, (msg)); } while (false)
#define REQUIRE_NEAR(a, b, tol) REQUIRE(::test::near((double)(a), (double)(b), (double)(tol)))
#define REQUIRE_VEC_NEAR(a, b, tol) REQUIRE(::test::nearVec((a), (b), (double)(tol)))

