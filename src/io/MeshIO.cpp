#include "io/MeshIO.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace bff_sdf::io {
namespace {

struct FaceToken {
    int v = -1;
    int vt = -1;
};

int parsePositiveObjIndex(const std::string& token, int count)
{
    if (token.empty()) return -1;
    const int raw = std::stoi(token);
    if (raw > 0) return raw - 1;
    if (raw < 0) return count + raw;
    return -1;
}

FaceToken parseFaceToken(const std::string& token, int vertexCount, int uvCount)
{
    FaceToken out;
    const size_t firstSlash = token.find('/');
    if (firstSlash == std::string::npos) {
        out.v = parsePositiveObjIndex(token, vertexCount);
        return out;
    }

    out.v = parsePositiveObjIndex(token.substr(0, firstSlash), vertexCount);
    const size_t secondSlash = token.find('/', firstSlash + 1);
    const std::string vt =
        secondSlash == std::string::npos
            ? token.substr(firstSlash + 1)
            : token.substr(firstSlash + 1, secondSlash - firstSlash - 1);
    if (!vt.empty()) out.vt = parsePositiveObjIndex(vt, uvCount);
    return out;
}

void setError(std::string* error, const std::string& message)
{
    if (error) *error = message;
}

} // namespace

bool readStl(const std::string& path, bff_sdf::core::Mesh& mesh, std::string* error)
{
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        setError(error, "unable to open STL: " + path);
        return false;
    }

    std::array<char, 80> header{};
    std::uint32_t triCount = 0;
    in.read(header.data(), static_cast<std::streamsize>(header.size()));
    in.read(reinterpret_cast<char*>(&triCount), sizeof(std::uint32_t));
    const auto fileSize = std::filesystem::file_size(path);
    const bool looksBinary = in && triCount > 0 && fileSize == 84ull + static_cast<std::uintmax_t>(triCount) * 50ull;

    std::vector<Eigen::Vector3d> vertices;
    std::vector<Eigen::Vector3i> faces;
    std::map<std::array<long long, 3>, int> welded;
    constexpr double scale = 1e9;

    auto addVertex = [&](const Eigen::Vector3d& p) {
        const std::array<long long, 3> key = {
            static_cast<long long>(std::llround(p.x() * scale)),
            static_cast<long long>(std::llround(p.y() * scale)),
            static_cast<long long>(std::llround(p.z() * scale))
        };
        auto it = welded.find(key);
        if (it != welded.end()) return it->second;
        const int id = static_cast<int>(vertices.size());
        welded[key] = id;
        vertices.push_back(p);
        return id;
    };

    if (looksBinary) {
        for (std::uint32_t t = 0; t < triCount; ++t) {
            float raw[12];
            std::uint16_t attr = 0;
            in.read(reinterpret_cast<char*>(raw), sizeof(raw));
            in.read(reinterpret_cast<char*>(&attr), sizeof(attr));
            if (!in) {
                setError(error, "truncated binary STL: " + path);
                return false;
            }
            Eigen::Vector3i f;
            for (int k = 0; k < 3; ++k) {
                const Eigen::Vector3d p(raw[3 + 3 * k], raw[4 + 3 * k], raw[5 + 3 * k]);
                f[k] = addVertex(p);
            }
            if (f[0] != f[1] && f[1] != f[2] && f[2] != f[0]) faces.push_back(f);
        }
    } else {
        in.close();
        std::ifstream text(path);
        std::string word;
        std::vector<int> tri;
        while (text >> word) {
            if (word == "vertex") {
                Eigen::Vector3d p;
                text >> p.x() >> p.y() >> p.z();
                tri.push_back(addVertex(p));
                if (tri.size() == 3) {
                    if (tri[0] != tri[1] && tri[1] != tri[2] && tri[2] != tri[0]) {
                        faces.emplace_back(tri[0], tri[1], tri[2]);
                    }
                    tri.clear();
                }
            }
        }
    }

    if (vertices.empty() || faces.empty()) {
        setError(error, "invalid or empty STL: " + path);
        return false;
    }

    mesh.V.resize(static_cast<int>(vertices.size()), 3);
    for (int i = 0; i < static_cast<int>(vertices.size()); ++i) mesh.V.row(i) = vertices[i].transpose();
    mesh.F.resize(static_cast<int>(faces.size()), 3);
    for (int i = 0; i < static_cast<int>(faces.size()); ++i) mesh.F.row(i) = faces[i].transpose();
    mesh.computeNormals();
    return true;
}

bool readObjWithUV(const std::string& path,
                   bff_sdf::core::Mesh& mesh,
                   Eigen::MatrixXd& UV,
                   Eigen::MatrixXi& Fuv,
                   std::string* error)
{
    std::ifstream in(path);
    if (!in) {
        setError(error, "unable to open OBJ: " + path);
        return false;
    }

    std::vector<Eigen::Vector3d> vertices;
    std::vector<Eigen::Vector2d> uvs;
    std::vector<Eigen::Vector3i> faces;
    std::vector<Eigen::Vector3i> uvFaces;
    std::string line;

    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ls(line);
        std::string tag;
        ls >> tag;
        if (tag == "v") {
            Eigen::Vector3d p;
            ls >> p.x() >> p.y() >> p.z();
            vertices.push_back(p);
        } else if (tag == "vt") {
            Eigen::Vector2d uv;
            ls >> uv.x() >> uv.y();
            uvs.push_back(uv);
        } else if (tag == "f") {
            std::vector<FaceToken> tokens;
            std::string token;
            while (ls >> token) tokens.push_back(parseFaceToken(token, (int)vertices.size(), (int)uvs.size()));
            if (tokens.size() < 3) continue;
            for (size_t k = 1; k + 1 < tokens.size(); ++k) {
                faces.emplace_back(tokens[0].v, tokens[k].v, tokens[k + 1].v);
                uvFaces.emplace_back(tokens[0].vt, tokens[k].vt, tokens[k + 1].vt);
            }
        }
    }

    mesh.V.resize((int)vertices.size(), 3);
    for (int i = 0; i < (int)vertices.size(); ++i) mesh.V.row(i) = vertices[i].transpose();
    mesh.F.resize((int)faces.size(), 3);
    for (int i = 0; i < (int)faces.size(); ++i) mesh.F.row(i) = faces[i].transpose();
    mesh.computeNormals();

    UV.resize((int)uvs.size(), 2);
    for (int i = 0; i < (int)uvs.size(); ++i) UV.row(i) = uvs[i].transpose();
    Fuv.resize((int)uvFaces.size(), 3);
    for (int i = 0; i < (int)uvFaces.size(); ++i) Fuv.row(i) = uvFaces[i].transpose();

    return true;
}

bool readObj(const std::string& path, bff_sdf::core::Mesh& mesh, std::string* error)
{
    Eigen::MatrixXd UV;
    Eigen::MatrixXi Fuv;
    return readObjWithUV(path, mesh, UV, Fuv, error);
}

bool writeObj(const std::string& path, const bff_sdf::core::Mesh& mesh, std::string* error)
{
    std::ofstream out(path);
    if (!out) {
        setError(error, "unable to write OBJ: " + path);
        return false;
    }
    for (int i = 0; i < mesh.V.rows(); ++i) {
        out << "v " << mesh.V(i, 0) << " " << mesh.V(i, 1) << " " << mesh.V(i, 2) << "\n";
    }
    for (int i = 0; i < mesh.F.rows(); ++i) {
        out << "f " << mesh.F(i, 0) + 1 << " " << mesh.F(i, 1) + 1 << " " << mesh.F(i, 2) + 1 << "\n";
    }
    return true;
}

bool writeObjWithUV(const std::string& path,
                    const bff_sdf::core::Mesh& mesh,
                    const Eigen::MatrixXd& UV,
                    const Eigen::MatrixXi& Fuv,
                    std::string* error)
{
    std::ofstream out(path);
    if (!out) {
        setError(error, "unable to write OBJ with UV: " + path);
        return false;
    }
    for (int i = 0; i < mesh.V.rows(); ++i) {
        out << "v " << mesh.V(i, 0) << " " << mesh.V(i, 1) << " " << mesh.V(i, 2) << "\n";
    }
    for (int i = 0; i < UV.rows(); ++i) {
        out << "vt " << UV(i, 0) << " " << UV(i, 1) << "\n";
    }
    for (int i = 0; i < mesh.F.rows(); ++i) {
        out << "f";
        for (int k = 0; k < 3; ++k) {
            const int vi = mesh.F(i, k) + 1;
            const int ti = Fuv.rows() == mesh.F.rows() ? Fuv(i, k) + 1 : vi;
            out << " " << vi << "/" << ti;
        }
        out << "\n";
    }
    return true;
}

} // namespace bff_sdf::io
